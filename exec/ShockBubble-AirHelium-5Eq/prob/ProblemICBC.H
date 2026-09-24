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
// Shock-helium bubble interaction (after Haas and Sturtevant, JFM 181, 1987). Helium
// (material 2) fills a circle of radius prob.r centered at the origin, on the lower slip
// wall, in air (material 1) at rest with pressure 1 and density prob.rho1. A shock of Mach
// number prob.Mach travels in +x through the air from x = -1.1 r. The interface is a tanh
// profile of width eta = PhaseField.Eta_Multiplier finest cells.
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
    // Interface: tanh profile of width eta about the bubble surface
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-10;
    const amrex::Real Psi = parm.r - std::sqrt(x*x + y*y + z*z);

    const amrex::Real Alpha1 = 1.0 - 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Gases at rest, and the post-shock state of the air behind the shock
    // ---------------------------------------
    const amrex::Real P_Unshocked = 1.0;
    const amrex::Real x_Shock     = -1.1 * parm.r;

    amrex::Real rho1 = parm.rho1;
    amrex::Real rho2 = parm.rho2;
    amrex::Real u    = 0.0;
    amrex::Real v    = 0.0;
    amrex::Real w    = 0.0;
    amrex::Real P    = P_Unshocked;

    if (x < x_Shock) {
        amrex::Real shock[5];
        StiffenedShockConditions(parm.Mach, 0.0, P_Unshocked, parm.rho1, parm.gamma_1, parm.pinf_1, shock);
        rho1 = shock[0];
        u    = shock[1];
        P    = shock[2];
    }

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1 * rho1, Alpha2 * rho2,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/bubble.csv with
//   V        helium volume, int alpha_2 dV (in 2D, the area of the upper half of the bubble)
//   x_c      x-centroid of the helium, int alpha_2 x dV / V
//   u_c      alpha_2-weighted mean x-velocity of the helium, int alpha_2 u dV / V
//   x_up     upstream interface on the axis, the smallest x where alpha_2 crosses 0.5
//   x_down   downstream interface on the axis, the largest x where alpha_2 crosses 0.5
// The integrals are over the composite grid; the axial crossings are taken on the finest
// level, in the row of cells next to y = 0 (and z = 0 in 3D).
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // One ghost cell, for the +x neighbor in the crossings
    const auto v = Measure(o, {Sum, Sum, Sum, Min, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a2  = 1.0 - U(i,j,k,INDEX_VolumeFraction1);
        const Real a2e = 1.0 - U(i+1,j,k,INDEX_VolumeFraction1);
        const Real u   = U(i,j,k,INDEX_MomentumX) / (U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2));

        // alpha_2 = 0.5 crossing toward the +x neighbor, on the axis
        bool axis = c.finest && (j == c.domain.smallEnd(1));
#if (AMREX_SPACEDIM == 3)
        axis = axis && (c.x[2] - 0.5 * c.dx[2] <= 0.0) && (c.x[2] + 0.5 * c.dx[2] > 0.0);
#endif
        Real x_up = big, x_down = -big;
        if (axis && (a2 - 0.5) * (a2e - 0.5) <= 0.0 && a2 != a2e) {
            x_up = x_down = c.x[0] + c.dx[0] * (0.5 - a2) / (a2e - a2);
        }

        return {a2 * c.dV, a2 * c.x[0] * c.dV, a2 * u * c.dV, x_up, x_down};
    });

    const Real V = v[0];
    WriteRow(o, "bubble.csv", {{"V", V}, {"x_c", (V > 0.0) ? v[1] / V : 0.0}, {"u_c", (V > 0.0) ? v[2] / V : 0.0},
                                {"x_up", v[3]}, {"x_down", v[4]}});
}

#endif  // IC_H_
