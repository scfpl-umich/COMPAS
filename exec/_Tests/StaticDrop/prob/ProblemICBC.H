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
// A drop of phase 2 (radius prob.R, center (prob.xc, prob.yc), density prob.rho_2) in phase 1
// (density prob.rho_1), at the velocity (prob.u0, prob.v0), 0 for the static drop. The volume
// fraction is a tanh profile of width eta about the signed distance, kept eps/2 away from 0 and
// 1. The pressure is prob.P0 outside and, with prob.Laplace_IC = 1, prob.P0 + sigma_12/R inside
// (2 sigma_12/R in 3D, a sphere), blended with the volume fraction of the drop; with
// prob.Laplace_IC = 0 it is uniform. The same file serves the four models; in the N-phase models
// the phases 3 to NPHASE are absent (volume fraction, partial density and energy exactly 0).
// =======================================

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
    amrex::Real r2 = (x - parm.xc)*(x - parm.xc) + (y - parm.yc)*(y - parm.yc);
#if (AMREX_SPACEDIM > 2)
    r2 += z*z;
#endif
    const amrex::Real psi = parm.R - std::sqrt(r2);
    const amrex::Real Alpha_Drop = 0.5*(1.0 + std::tanh( psi/(std::sqrt(2.0)*eta) )*(1.0 - eps) );

#if (SURFACE_TENSION == true)
    const amrex::Real sigma = parm.Physics_Parm.SurfaceTension_Parm.sigma_12;
#else
    const amrex::Real sigma = 0.0;
#endif
    const amrex::Real P = parm.P0 + parm.Laplace_IC*(AMREX_SPACEDIM - 1)*sigma/parm.R*Alpha_Drop;

    const amrex::Real u = parm.u0;
    const amrex::Real v = parm.v0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

    amrex::Array<amrex::Real, NSTATE> Primitive;
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    const amrex::Real Alpha1 = 1.0 - Alpha_Drop;
    const amrex::Real Alpha2 = Alpha_Drop;
#if (PHYSICS == FIVEEQS)
    State2Primitive(Primitive,
                    Alpha1*parm.rho_1, Alpha2*parm.rho_2,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
#else
    State2Primitive(Primitive,
                    Alpha1*parm.rho_1, Alpha2*parm.rho_2,
       AMREX_D_DECL(u, v, w),
                    Alpha1*P, Alpha2*P,
                    Alpha1);
#endif
#else
    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha[i_Phase] = 0.0;
    }
    Alpha[1] = Alpha_Drop;
    Alpha[0] = 1.0 - Alpha_Drop;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        const amrex::Real rho_k = (i_Phase == 0) ? parm.rho_1 : ((i_Phase == 1) ? parm.rho_2 : 0.0);
        Alpha_rho[i_Phase] = Alpha[i_Phase]*rho_k;
        Alpha_P  [i_Phase] = Alpha[i_Phase]*P;
    }
#if (PHYSICS == FIVEEQS_NPHASE)
    State2Primitive(Primitive,
                    Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha);
#else
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    P);
#endif
#endif
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends one
// line to <plot directory>/drop.csv with
//   P_in, P_out        mean mixture pressure of the cells with alpha_2 > 0.999 (inside the
//                      drop) and alpha_2 < 0.001 (outside)
//   dP, dP_exact       P_in - P_out and the Laplace jump (AMREX_SPACEDIM - 1) sigma_12/R
//   u_max              max |u| (spurious currents; Ca = mu |u|_max/sigma)
//   mom_x, mom_y, mom_z  int rho u dV, the momentum (exactly conserved in a periodic box)
//   mom_abs            int |rho u| dV, the scale of the momentum error
//   mass, energy       int rho dV and the total energy int rho E dV
//   kinetic            int rho |u|^2/2 dV
//   V_2                int alpha_2 dV, the drop volume (area in 2D)
//   absent             max over the phases 3 to NPHASE of |alpha_k|, |alpha_k rho_k| and
//                      their energy, which must stay exactly 0 (0 for the two-phase models)
// The integrals and the extrema are over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto r = Measure(o, {Sum, Sum, Sum, Sum, Max, Sum, Sum, Sum, Sum, Sum, Sum, Sum, Sum, Max}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<14>
    {
        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Conservative2Primitive(Prim, Cons, parm);

        Real rho = 0.0;
        Real Alpha_Drop;
        Real P;
        Real E;
        Real absent = 0.0;
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
        rho = Cons[INDEX_Mass1] + Cons[INDEX_Mass2];
        Alpha_Drop = 1.0 - Cons[INDEX_VolumeFraction1];
#if (PHYSICS == FIVEEQS)
        P = Prim[INDEX_Energy];
        E = Cons[INDEX_Energy];
#else
        P = Prim[INDEX_Energy1] + Prim[INDEX_Energy2];
        E = Cons[INDEX_Energy1] + Cons[INDEX_Energy2];
#endif
#else
        for (int n = 0; n < NPHASE; ++n) { rho += Cons[INDEX_Mass1+n]; }
        Alpha_Drop = Cons[INDEX_VolumeFraction1+1];
#if (PHYSICS == FIVEEQS_NPHASE)
        P = Prim[INDEX_Energy];
        E = Cons[INDEX_Energy];
        for (int n = 2; n < NPHASE; ++n) {
            absent = std::max({absent, std::abs(Cons[INDEX_VolumeFraction1+n]), std::abs(Cons[INDEX_Mass1+n])});
        }
#else
        P = Prim[INDEX_Total_Energy];
        E = Cons[INDEX_Total_Energy];
        for (int n = 2; n < NPHASE; ++n) {
            absent = std::max({absent, std::abs(Cons[INDEX_VolumeFraction1+n]),
                               std::abs(Cons[INDEX_Mass1+n]), std::abs(Cons[INDEX_Energy1+n])});
        }
#endif
#endif
        amrex::Array<Real, 3> Mom = {0.0, 0.0, 0.0};
        Real Vel2 = 0.0;
        Real Mom_abs = 0.0;
        for (int b = 0; b < AMREX_SPACEDIM; ++b) {
            Mom[b] = Cons[INDEX_MomentumX+b];
            Vel2 += (Mom[b]/rho)*(Mom[b]/rho);
            Mom_abs += std::abs(Mom[b]);
        }
        const bool Inside  = Alpha_Drop > 0.999;
        const bool Outside = Alpha_Drop < 0.001;

        return {Inside ? P*c.dV : 0.0, Inside ? c.dV : 0.0,
                Outside ? P*c.dV : 0.0, Outside ? c.dV : 0.0,
                std::sqrt(Vel2),
                Mom[0]*c.dV, Mom[1]*c.dV, Mom[2]*c.dV, Mom_abs*c.dV,
                rho*c.dV, E*c.dV, 0.5*rho*Vel2*c.dV, Alpha_Drop*c.dV, absent};
    });

    const Real P_in  = (r[1] > 0.0) ? r[0]/r[1] : 0.0;
    const Real P_out = (r[3] > 0.0) ? r[2]/r[3] : 0.0;
#if (SURFACE_TENSION == true)
    const Real sigma = HostParm(o).Physics_Parm.SurfaceTension_Parm.sigma_12;
#else
    const Real sigma = 0.0;
#endif
    const Real dP_exact = (AMREX_SPACEDIM - 1)*sigma/HostParm(o).R;

    WriteRow(o, "drop.csv", {{"P_in", P_in}, {"P_out", P_out}, {"dP", P_in - P_out}, {"dP_exact", dP_exact},
                             {"u_max", r[4]}, {"mom_x", r[5]}, {"mom_y", r[6]}, {"mom_z", r[7]},
                             {"mom_abs", r[8]}, {"mass", r[9]}, {"energy", r[10]}, {"kinetic", r[11]},
                             {"V_2", r[12]}, {"absent", r[13]}});
}

#endif  // IC_H_
