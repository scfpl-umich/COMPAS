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
// Planar water hammer. Water (material 1, stiffened gas) with pressure prob.P0 and density
// prob.rho0 fills the strip and moves in +x at prob.u_in toward the slip wall at the right
// end. Material 2 is present only as a trace volume fraction prob.eps, with density
// prob.rho2, in pressure equilibrium with the water. The state is uniform, so the wall
// starts the reflected shock at t = 0; its exact speed and the water-hammer pressure
// behind it are computed in Parm::DynamicInit. SI units.
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
    amrex::ignore_unused(x, y, z);
    amrex::ignore_unused(AMREX_D_DECL(dx, dy, dz));
    amrex::ignore_unused(AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    // ---------------------------------------
    // Incoming water with the trace of material 2
    // ---------------------------------------
    const amrex::Real Alpha1 = 1.0 - parm.eps;
    const amrex::Real Alpha2 = parm.eps;

    const amrex::Real u = parm.u_in;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(u, v, w);

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*parm.rho0, Alpha2*parm.rho2,
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
// one line to <plot directory>/waterhammer.csv with
//   x_shock   position of the reflected shock, the smallest x where the pressure crosses
//             the mean of prob.P0 and the exact water-hammer pressure, found on the finest
//             level (the wall position before the shock has formed)
//   P_wall    pressure in the cells next to the wall, averaged over the strip
//   u_wall    velocity along x in the cells next to the wall, averaged over the strip
// The wall cells are taken on the finest level that covers them.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // One ghost cell, for the crossing with the +x neighbor
    const auto v = Measure(o, {Min, Sum, Sum, Sum}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<4>
    {
        amrex::Array<Real, NSTATE> Q, Qe;
        for (int n = 0; n < NSTATE; ++n) { Q[n] = U(c.i,c.j,c.k,n); Qe[n] = U(c.i+1,c.j,c.k,n); }
        Real P, Pe;
        Conservative2Pressure(P,  Q,  parm);
        Conservative2Pressure(Pe, Qe, parm);

        // Crossing of the mid pressure toward the +x neighbor
        const Real P_mid = 0.5 * (parm.P0 + parm.P_wh);
        const bool cross = c.finest && (P - P_mid) * (Pe - P_mid) <= 0.0 && P != Pe;
        const Real x_s   = cross ? c.x[0] + c.dx[0] * (P_mid - P) / (Pe - P) : big;

        // Cells next to the wall, on the finest level that covers them
        const bool wall = (c.i == c.domain.bigEnd(0));
        const Real rho  = Q[INDEX_Mass1] + Q[INDEX_Mass2];
        return {x_s, wall ? c.dV : 0.0, wall ? P * c.dV : 0.0, wall ? Q[INDEX_MomentumX] / rho * c.dV : 0.0};
    });

    const Real x_wall  = FinestGeom(o).ProbHi(0);
    const Real x_shock = (v[0] < x_wall) ? v[0] : x_wall;
    const Real A_wall  = v[1];
    WriteRow(o, "waterhammer.csv", {{"x_shock", x_shock},
                                    {"P_wall", (A_wall > 0.0) ? v[2] / A_wall : 0.0},
                                    {"u_wall", (A_wall > 0.0) ? v[3] / A_wall : 0.0}});
}

#endif  // IC_H_
