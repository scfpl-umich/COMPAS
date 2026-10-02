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
// Viscous shock tube (Daru and Tenaud, Comput. Fluids 38, 2009). Gas at rest with a
// diaphragm at x = 0.5: rho = 120, P = 120/gamma on the left and rho = 1.2, P = 1.2/gamma
// on the right. One gas: material 1 fills the box, with a trace eps of material 2, the same
// gas with the same transport properties, so the result can be compared with
// ViscousShockTube-5Eq.
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
    amrex::ignore_unused(y, z, AMREX_D_DECL(dx, dy, dz));

    // ---------------------------------------
    // Left and right states
    // ---------------------------------------
    const amrex::Real gamma = parm.gamma_1;
    const bool left = (x < 0.5);
    const amrex::Real rho = left ? 120.0 : 1.2;
    const amrex::Real p   = rho / gamma;
    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    const amrex::Real eps = 1.0e-6; // Trace of the second phase
    const amrex::Real Alpha1 = 1.0 - eps;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*rho, (1.0 - Alpha1)*rho,
       AMREX_D_DECL(u, v, w),
                    Alpha1*p, (1.0 - Alpha1)*p,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/box.csv with the integrals over the composite grid of
//   mass       density
//   energy     total energy, alpha_1 rho_1 E_1 + alpha_2 rho_2 E_2
//   kinetic    kinetic energy
// The box is closed by adiabatic walls, so mass and total energy are conserved while the
// kinetic energy changes through compression and viscous dissipation.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto v = Measure<3>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<3>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real rho = U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2);
        const Real m2  = AMREX_D_TERM(U(i,j,k,INDEX_MomentumX) * U(i,j,k,INDEX_MomentumX),
                                    + U(i,j,k,INDEX_MomentumY) * U(i,j,k,INDEX_MomentumY),
                                    + U(i,j,k,INDEX_MomentumZ) * U(i,j,k,INDEX_MomentumZ));
        return {rho * c.dV, (U(i,j,k,INDEX_Energy1) + U(i,j,k,INDEX_Energy2)) * c.dV, 0.5 * m2 / rho * c.dV};
    });
    WriteRow(o, "box.csv", {{"mass", v[0]}, {"energy", v[1]}, {"kinetic", v[2]}});
}

#endif  // IC_H_
