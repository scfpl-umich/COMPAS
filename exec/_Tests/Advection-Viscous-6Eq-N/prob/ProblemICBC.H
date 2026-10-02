// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>

static_assert(NPHASE >= 2, "This case needs at least two phases.");


// =======================================
// Initial condition
//
// A liquid drop of phase 2 (radius prob.R, center (prob.xc, prob.yc), density prob.rho_2) in a
// gas, phase 1 (density prob.rho_1), at the uniform pressure prob.P0 and velocity
// (prob.u0, prob.v0). The volume fraction is a tanh profile of width eta about the signed
// distance, kept eps/2 away from 0 and 1. Phases 3 to NPHASE are absent: their volume
// fraction, partial density and internal energy are exactly 0, and they must stay so. Each
// present phase has a uniform density, so its temperature is uniform too, and the two
// temperatures differ. With viscosity and heat conduction, the exact solution is still a
// translation: the velocity is uniform, so the stress vanishes, and each phase temperature is
// uniform, so no phase conducts heat (Abgrall, JCP 125, 1996).
// =======================================

// Temperature of phase ID_Phase (0-based) at the pressure P0 and the density rho, from its EOS
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real DropTemperature(int ID_Phase, amrex::Real rho, Parm const& parm)
{
    const amrex::Real A  = parm.A[ID_Phase];
    const amrex::Real B  = parm.B[ID_Phase];
    const amrex::Real cv = parm.cv[ID_Phase];
    return (A*parm.P0 + A*B/(1.0 + A))/(cv*rho);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
{
    // Prevent unused-variable warnings
    sharp = false;
    amrex::ignore_unused(z, AMREX_D_DECL(dx, dy, dz));

    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;
    const amrex::Real psi = parm.R - std::sqrt((x - parm.xc)*(x - parm.xc) + (y - parm.yc)*(y - parm.yc));

    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha[i_Phase] = 0.0;
    }
    Alpha[1] = 0.5*(1.0 + std::tanh( psi/(std::sqrt(2.0)*eta) )*(1.0 - eps) );
    Alpha[0] = 1.0 - Alpha[1];
    const amrex::Real P = parm.P0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        const amrex::Real rho_k = (i_Phase == 0) ? parm.rho_1 : ((i_Phase == 1) ? parm.rho_2 : 0.0);
        Alpha_rho[i_Phase] = Alpha[i_Phase]*rho_k;
        Alpha_P  [i_Phase] = Alpha[i_Phase]*P;
    }

    const amrex::Real u = parm.u0;
    const amrex::Real v = parm.v0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    P);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/advection.csv with
//   V_2                   volume (area in 2D) of the drop, int alpha_2 dV
//   alpha_min, alpha_max  extrema of alpha_1
//   dP                    max |P - P0|, with the mixture pressure P = sum_k alpha_k P_k
//   du                    max |u - u0| and |v - v0|
//   dT_1, dT_2            max |T_k - T_k0|/T_k0 over the cells with alpha_k > 0.01, with T_k0 the
//                         initial temperature of phase k
//   absent                max over the phases 3 to NPHASE of |alpha_k|, |alpha_k rho_k| and
//                         |alpha_k rho_k e_k|, which must stay exactly 0
//   mass_1, mass_2        int alpha_k rho_k dV
//   energy                int rho E dV
// The integrals and the extrema are over the composite grid. Exact advection keeps the uniform
// state, so dP, du, dT_1 and dT_2 measure the errors.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto r = Measure(o, {Sum, Min, Max, Max, Max, Max, Max, Max, Sum, Sum, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<11>
    {
        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Conservative2Primitive(Prim, Cons, parm);
        amrex::Array<Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
        Real P;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(Alpha, Alpha_rho, Alpha_P, AMREX_D_DECL(u, v, w), P, Prim);

        amrex::Array<Real, NPHASE> T_k;
        amrex::Array<bool, NPHASE> Present;
        Conservative2PhaseTemperature(T_k, Present, Cons, parm);
        const Real T10 = DropTemperature(0, parm.rho_1, parm);
        const Real T20 = DropTemperature(1, parm.rho_2, parm);
        const Real dT1 = (Alpha[0] > 0.01) ? std::abs(T_k[0] - T10)/T10 : 0.0;
        const Real dT2 = (Alpha[1] > 0.01) ? std::abs(T_k[1] - T20)/T20 : 0.0;

        Real absent = 0.0;
        for (int n = 2; n < NPHASE; ++n) {
            absent = std::max({absent, std::abs(Cons[INDEX_VolumeFraction1+n]),
                               std::abs(Cons[INDEX_Mass1+n]), std::abs(Cons[INDEX_Energy1+n])});
        }

        return {Alpha[1] * c.dV, Alpha[0], Alpha[0], std::abs(P - parm.P0),
                std::max(std::abs(u - parm.u0), std::abs(v - parm.v0)), dT1, dT2, absent,
                Alpha_rho[0] * c.dV, Alpha_rho[1] * c.dV, Cons[INDEX_Total_Energy] * c.dV};
    });

    WriteRow(o, "advection.csv", {{"V_2", r[0]}, {"alpha_min", r[1]}, {"alpha_max", r[2]},
                                  {"dP", r[3]}, {"du", r[4]}, {"dT_1", r[5]}, {"dT_2", r[6]},
                                  {"absent", r[7]}, {"mass_1", r[8]}, {"mass_2", r[9]},
                                  {"energy", r[10]}});
}

#endif  // IC_H_
