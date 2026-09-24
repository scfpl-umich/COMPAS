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
// Sod shock tube (Sod, J. Comput. Phys. 27, 1978) along direction prob.direction
// (0 = x, 1 = y, 2 = z), with the diaphragm at 0.5 and the gas at rest. Material 1
// fills the tube, with a trace volume fraction eps of material 2.
//
// These are the defaults. The optional parameters prob.rho_L, prob.u_L, prob.P_L,
// prob.rho_R, prob.u_R, prob.P_R, prob.x_diaphragm, prob.mat_L and prob.mat_R (see
// Parm.H) set the two states, the diaphragm position and the material (1 or 2) on each
// side, for example for a two-material shock tube. The other material is the trace eps.
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
    amrex::ignore_unused(AMREX_D_DECL(dx, dy, dz));

    // ---------------------------------------
    // Left and right states
    // ---------------------------------------
    const amrex::Real xi   = (parm.direction == 0) ? x : ((parm.direction == 1) ? y : z);
    const bool        left = (xi < parm.x_diaphragm);
    const amrex::Real rho  = left ? parm.rho_L : parm.rho_R;
    const amrex::Real P    = left ? parm.P_L   : parm.P_R;
    const amrex::Real uxi  = left ? parm.u_L   : parm.u_R;   // Velocity along the tube
    const int         mat  = left ? parm.mat_L : parm.mat_R; // Material that fills this side

    const amrex::Real u = (parm.direction == 0) ? uxi : 0.0;
    const amrex::Real v = (parm.direction == 1) ? uxi : 0.0;
    const amrex::Real w = (parm.direction == 2) ? uxi : 0.0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    const amrex::Real eps = 1.0e-6; // Trace of the second phase
    const amrex::Real Alpha1 = (mat == 1) ? 1.0 - eps : eps;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*rho, (1.0 - Alpha1)*rho,
       AMREX_D_DECL(u, v, w),
                    Alpha1*P, (1.0 - Alpha1)*P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/shock.csv with
//   x_shock    shock position along the tube, the furthest point where the density crosses
//              1.5 times the density ahead of the shock (prob.rho_R), found on the finest level
//   mass       total mass
//   momentum   total momentum along the tube
// The integrals are over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const int  nmat    = (NPHASE == 1) ? 2 : NPHASE;   // number of partial densities

    // One ghost cell, for the crossing with the next cell along the tube
    const auto v = Measure(o, {Max, Sum, Sum}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<3>
    {
        const int  d       = parm.direction;
        const Real rho_thr = 1.5 * parm.rho_R;
        const int ie = c.i + (d == 0), je = c.j + (d == 1), ke = c.k + (d == 2);
        Real rho = 0.0, rho_e = 0.0;
        for (int n = 0; n < nmat; ++n) { rho += U(c.i,c.j,c.k,INDEX_Mass1+n); rho_e += U(ie,je,ke,INDEX_Mass1+n); }

        const bool cross = c.finest && (rho - rho_thr) * (rho_e - rho_thr) <= 0.0 && rho != rho_e;
        return {cross ? c.x[d] + c.dx[d] * (rho_thr - rho) / (rho_e - rho) : -1.0e30,
                rho * c.dV, U(c.i,c.j,c.k,INDEX_MomentumX+d) * c.dV};
    });
    WriteRow(o, "shock.csv", {{"x_shock", v[0]}, {"mass", v[1]}, {"momentum", v[2]}});
}

#endif  // IC_H_
