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
// Decay of a shear wave u = U sin(k (x_a - a_lo)), v = w = 0, with x_a the coordinate along the
// axis prob.axis (1: y, or 2: z in 3D), a_lo its lower bound, k = 2 pi prob.mode/L_a and L_a the
// length of the box along that axis, in a uniform mixture of the two phases (volume fraction
// prob.alpha_1 of phase 1, densities prob.rho_1 and prob.rho_2) at the uniform pressure prob.P0.
// With the mixture viscosity mu = alpha_1 mu_1 + alpha_2 mu_2 and the density rho, the
// amplitude decays as exp(-nu k^2 t), nu = mu/rho, as long as the Mach number U/c is small.
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
    amrex::ignore_unused(x, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    const amrex::Real k   = 2.0*PI*parm.mode/parm.L_a;
    const amrex::Real x_a = (parm.axis == 2) ? z : y;
    const amrex::Real u   = parm.U*std::sin(k*(x_a - parm.a_lo));
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

    const amrex::Real Alpha1 = parm.alpha_1;
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*parm.rho_1, Alpha2*parm.rho_2,
       AMREX_D_DECL(u, v, w),
                    Alpha1*parm.P0, Alpha2*parm.P0,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/shear.csv with the integrals over the composite grid of
//   mass_1, mass_2     partial densities alpha_k rho_k
//   mom_x, mom_y       momentum
//   energy             total energy, alpha_1 rho_1 E_1 + alpha_2 rho_2 E_2
//   kinetic            kinetic energy rho |u|^2/2
//   ie_1, ie_2         phase internal energies alpha_k rho_k e_k = alpha_k rho_k (E_k - |u|^2/2)
// and the amplitude of the wave, amp = 2/V int u sin(k (x_a - a_lo)) dV, with V the volume of
// the box. The box is periodic, so mass, momentum and energy are conserved, while the kinetic
// energy goes into the internal energies, in the ratio alpha_1 mu_1 : alpha_2 mu_2 without relaxation.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto v = Measure<9>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<9>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1r1 = U(i,j,k,INDEX_Mass1);
        const Real a2r2 = U(i,j,k,INDEX_Mass2);
        const Real rho  = a1r1 + a2r2;
        AMREX_D_TERM(const Real ru = U(i,j,k,INDEX_MomentumX);,
                     const Real rv = U(i,j,k,INDEX_MomentumY);,
                     const Real rw = U(i,j,k,INDEX_MomentumZ);)
        const Real m2 = AMREX_D_TERM(ru*ru, + rv*rv, + rw*rw);
        const Real v2 = m2/(rho*rho);
        const Real E1 = U(i,j,k,INDEX_Energy1);
        const Real E2 = U(i,j,k,INDEX_Energy2);
        const Real kw = 2.0*PI*parm.mode/parm.L_a;
#if (AMREX_SPACEDIM > 1)
        const Real my = rv;
#else
        const Real my = 0.0;
#endif
        return {a1r1*c.dV, a2r2*c.dV, ru*c.dV, my*c.dV, (E1 + E2)*c.dV, 0.5*m2/rho*c.dV,
                (E1 - 0.5*a1r1*v2)*c.dV, (E2 - 0.5*a2r2*v2)*c.dV,
                2.0/parm.Volume*ru/rho*std::sin(kw*(c.x[AMREX_SPACEDIM > 1 ? parm.axis : 0] - parm.a_lo))*c.dV};
    });
    WriteRow(o, "shear.csv", {{"mass_1", v[0]}, {"mass_2", v[1]}, {"mom_x", v[2]}, {"mom_y", v[3]},
                              {"energy", v[4]}, {"kinetic", v[5]}, {"ie_1", v[6]}, {"ie_2", v[7]},
                              {"amp", v[8]}});
}

#endif  // IC_H_
