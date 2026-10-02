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


// =======================================
// Initial condition
//
// Tests of the Noble-Abel stiffened-gas EOS (prob.EOS.EOS = 3; Le Metayer and Saurel, Phys.
// Fluids 28, 046102, 2016), in SI units, selected by prob.problem:
//
// 0  Shock tube along x: material prob.mat_L at (prob.rho_L, prob.u_L, prob.P_L) for
//    x < prob.x_diaphragm and material prob.mat_R at (prob.rho_R, prob.u_R, prob.P_R) to the
//    right. The other material is a trace (alpha = 1e-6) with the same density and pressure.
//    The base inputs are air (ideal gas) at 1e9 Pa against liquid water (NASG) at 1e5 Pa.
// 1  Drop: a disk of material 1 of radius prob.R_drop centered at the origin, in material 2,
//    as a tanh profile of width eta = PhaseField.Eta_Multiplier finest cells, kept eps/2
//    (eps = 1e-12) away from 0 and 1. Both materials are at prob.P_drop and prob.T_drop and
//    move with (prob.u_drop, prob.u_drop). The phase densities follow from the EOS,
//    rho_k = 1/(b_k + (gamma_k - 1) cv_k T/(P + pinf_k)), so they are uniform, in the traces
//    too, and exact advection keeps the pressure and the velocity uniform (Abgrall, JCP 125, 1996).
// =======================================

// Phase density at pressure P and temperature T (NASG; the stiffened gas for b = 0)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real PhaseDensity_PT (amrex::Real P, amrex::Real T, amrex::Real gamma, amrex::Real pinf,
                             amrex::Real b, amrex::Real cv)
{
    return 1.0/(b + (gamma - 1.0)*cv*T/(P + pinf));
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
    amrex::ignore_unused(y, z, AMREX_D_DECL(dx, dy, dz));

    amrex::Array<amrex::Real, NSTATE> Primitive;
    if (parm.problem == 1){
        // ---------------------------------------
        // Drop
        // ---------------------------------------
        const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
        const amrex::Real eps = 1.0e-12;
        const amrex::Real psi = parm.R_drop - std::sqrt(AMREX_D_TERM(x*x, + y*y, + z*z));
        const amrex::Real Alpha1 = 0.5*(1.0 + std::tanh( psi/(std::sqrt(2.0)*eta) )*(1.0 - eps));
        const amrex::Real Alpha2 = 1.0 - Alpha1;
        const amrex::Real P  = parm.P_drop;
        const amrex::Real U0 = parm.u_drop;
        amrex::Array<amrex::Real, NPHASE> Alpha = {Alpha1, Alpha2};
        amrex::Array<amrex::Real, NPHASE> Alpha_rho;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Alpha_rho[i_Phase] = Alpha[i_Phase]*PhaseDensity_PT(P, parm.T_drop, parm.gamma[i_Phase], parm.pinf[i_Phase],
                                                                parm.C[i_Phase], parm.cv[i_Phase]);
        }
        amrex::Array<amrex::Real, NPHASE> Alpha_P = {Alpha1*P, Alpha2*P};
        State2Primitive(Primitive,
                        Alpha, Alpha_rho, Alpha_P,
           AMREX_D_DECL(U0, U0, 0.0),
                        P);
    }
    else {
        // ---------------------------------------
        // Shock tube
        // ---------------------------------------
        const bool        left = (x < parm.x_diaphragm);
        const amrex::Real rho  = left ? parm.rho_L : parm.rho_R;
        const amrex::Real P    = left ? parm.P_L   : parm.P_R;
        const amrex::Real u    = left ? parm.u_L   : parm.u_R;
        const int         mat  = left ? parm.mat_L : parm.mat_R;
        const amrex::Real eps = 1.0e-6; // Trace of the other material
        amrex::Array<amrex::Real, NPHASE> Alpha = (mat == 1) ? amrex::Array<amrex::Real, NPHASE>{1.0 - eps, eps}
                                                             : amrex::Array<amrex::Real, NPHASE>{eps, 1.0 - eps};
        amrex::Array<amrex::Real, NPHASE> Alpha_rho;
        amrex::Array<amrex::Real, NPHASE> Alpha_P;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Alpha_rho[i_Phase] = Alpha[i_Phase]*rho;
            Alpha_P  [i_Phase] = Alpha[i_Phase]*P;
        }
        State2Primitive(Primitive,
                        Alpha, Alpha_rho, Alpha_P,
           AMREX_D_DECL(u, 0.0, 0.0),
                        P);
    }
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends one
// line to a file of the plot directory; the integrals are over the composite grid.
//
// Shock tube, nasg.csv:
//   mass_1, mass_2   int alpha_k rho_k dV, conserved
//   momentum         int rho u dV, which changes at the rate (P(0) - P(1)) L_y through the
//                    outflow boundaries while the waves are inside
//   energy           int rho E dV, conserved while the waves are inside
//   xi_min           smallest covolume factor 1 - b_k rho_k of the phases with alpha_k > 1e-6
//   P_min, P_max     extrema of the mixture pressure sum_k alpha_k p_k
//   x_contact        position of the alpha_1 = 1/2 crossing on the finest level
// Drop, drop.csv:
//   V_1                   int alpha_1 dV
//   mass_1, mass_2        int alpha_k rho_k dV
//   momentum_x, _y        int rho u dV, int rho v dV
//   energy                as above
//   alpha_min, alpha_max  extrema of alpha_1
//   xi_min                as above
//   dP                    max |P - P_drop|, with P the mixture pressure sum_k alpha_k p_k
//   du                    max |u - u_drop| and |v - u_drop|
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    if (HostParm(o).problem == 1){
        const auto r = Measure(o, {Sum, Sum, Sum, Sum, Sum, Sum, Min, Max, Min, Max, Max}, 0,
        [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<11>
        {
            amrex::Array<Real, NSTATE> Cons, Prim;
            for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
            Conservative2Primitive(Prim, Cons, parm);
            amrex::Array<Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
            Real P;
            AMREX_D_TERM(Real u;, Real v;, Real w;)
            Primitive2State(Alpha, Alpha_rho, Alpha_P, AMREX_D_DECL(u, v, w), P, Prim);
            amrex::ignore_unused(AMREX_D_DECL(u, v, w));
            const Real a1r1 = Alpha_rho[0], a2r2 = Alpha_rho[1], a1 = Alpha[0];
            Real xi = 1.0;
            for (int n = 0; n < NPHASE; ++n) {
                if (Alpha[n] > 1.0e-6) { xi = amrex::min(xi, 1.0 - parm.C[n]*Alpha_rho[n]/Alpha[n]); }
            }
            const Real du = std::max(std::abs(u - parm.u_drop), std::abs(v - parm.u_drop));
            return {a1 * c.dV, a1r1 * c.dV, a2r2 * c.dV, U(c.i,c.j,c.k,INDEX_MomentumX) * c.dV,
                    U(c.i,c.j,c.k,INDEX_MomentumY) * c.dV, U(c.i,c.j,c.k,INDEX_Total_Energy) * c.dV,
                    a1, a1, xi, std::abs(P - parm.P_drop), du};
        });
        WriteRow(o, "drop.csv", {{"V_1", r[0]}, {"mass_1", r[1]}, {"mass_2", r[2]}, {"momentum_x", r[3]},
                                 {"momentum_y", r[4]}, {"energy", r[5]}, {"alpha_min", r[6]},
                                 {"alpha_max", r[7]}, {"xi_min", r[8]}, {"dP", r[9]}, {"du", r[10]}});
        return;
    }

    // One ghost cell, for the crossing with the next cell along x
    const auto r = Measure(o, {Sum, Sum, Sum, Sum, Min, Min, Max, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<8>
    {
        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Conservative2Primitive(Prim, Cons, parm);
        amrex::Array<Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
        Real P;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(Alpha, Alpha_rho, Alpha_P, AMREX_D_DECL(u, v, w), P, Prim);
        amrex::ignore_unused(AMREX_D_DECL(u, v, w));
        const Real a1r1 = Alpha_rho[0], a2r2 = Alpha_rho[1], a1 = Alpha[0];
        Real xi = 1.0;
        for (int n = 0; n < NPHASE; ++n) {
            if (Alpha[n] > 1.0e-6) { xi = amrex::min(xi, 1.0 - parm.C[n]*Alpha_rho[n]/Alpha[n]); }
        }

        const Real a1_e  = U(c.i+1,c.j,c.k,INDEX_VolumeFraction1);
        const bool cross = c.finest && (a1 - 0.5) * (a1_e - 0.5) <= 0.0 && a1 != a1_e;

        return {a1r1 * c.dV, a2r2 * c.dV, U(c.i,c.j,c.k,INDEX_MomentumX) * c.dV,
                U(c.i,c.j,c.k,INDEX_Total_Energy) * c.dV, xi, P, P,
                cross ? c.x[0] + c.dx[0] * (0.5 - a1) / (a1_e - a1) : -1.0e30};
    });
    WriteRow(o, "nasg.csv", {{"mass_1", r[0]}, {"mass_2", r[1]}, {"momentum", r[2]}, {"energy", r[3]},
                             {"xi_min", r[4]}, {"P_min", r[5]}, {"P_max", r[6]}, {"x_contact", r[7]}});
}

#endif  // IC_H_
