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
// Four shapes of material 2 in material 1: Zalesak's slotted disk (Zalesak, JCP 31, 1979),
// a square, a disk and an equilateral triangle. Each shape is a tanh profile of width
// eta = Eta_Multiplier * dx_FinestLevel about its signed-distance function, and each has
// its own density (2, 3, 4 and 5; material 1 has density 1). Pressure is uniform and the
// velocity is (1, 1), so the exact solution is a translation of the initial condition.
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

    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;

    // Smoothed indicator of a shape from its signed distance (positive inside)
    auto indicator = [&](amrex::Real psi)
    {
        return 0.5*(1.0 + std::tanh(psi/(std::sqrt(2.0)*eta))*(1.0 - eps));
    };

    // ---------------------------------------
    // Shape 1: Zalesak's slotted disk, center (0.4, -0.4), radius 0.4, slot of width 0.125
    // open at the bottom and reaching 0.25 above the center
    // ---------------------------------------
    amrex::Real X = x - 0.4;
    amrex::Real Y = y + 0.4;
    const amrex::Real r_disk   = 0.4;
    const amrex::Real slot_w   = 0.125;
    const amrex::Real slot_d   = 0.65;
    const amrex::Real psi_disk = r_disk - std::sqrt(X*X + Y*Y);
    const amrex::Real psi_slot = std::min(0.5*slot_w - std::abs(X), (slot_d - r_disk) - Y);
    const amrex::Real H1 = indicator(std::min(psi_disk, -psi_slot));

    // ---------------------------------------
    // Shape 2: square of side 0.4, center (0.4, 0.5)
    // ---------------------------------------
    const amrex::Real H2 = indicator(0.2 - std::max(std::abs(x - 0.4), std::abs(y - 0.5)));

    // ---------------------------------------
    // Shape 3: disk of radius 0.2, center (-0.5, -0.5)
    // ---------------------------------------
    const amrex::Real H3 = indicator(0.2 - std::sqrt((x + 0.5)*(x + 0.5) + (y + 0.5)*(y + 0.5) + z*z));

    // ---------------------------------------
    // Shape 4: equilateral triangle of side 0.5 pointing in -y, with its upper edge at
    // y = 0.5 + h/3 and its lower vertex at (-0.5, 0.5 - h/3), where h = (sqrt(3)/2) 0.75
    // ---------------------------------------
    X = x + 0.5;
    Y = y - 0.5;
    const amrex::Real h  = std::sqrt(3.0) * 0.5 * 0.75;
    const amrex::Real d1 =  (h/3.0) - Y;
    const amrex::Real d2 =  (std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    const amrex::Real d3 = -(std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    const amrex::Real H4 = indicator(std::min(d1, std::min(d2, d3)));

    // ---------------------------------------
    // Material 2 fills the shapes, with the density of the shape it belongs to
    // ---------------------------------------
    const amrex::Real Alpha2 = H1 + H2 + H3 + H4;
    const amrex::Real Alpha1 = 1.0 - Alpha2;

    auto inside = [](amrex::Real H) { return (H > 1.0e-5) ? 1.0 : 0.0; };
    const amrex::Real rho1 = 1.0;
    const amrex::Real rho2 = inside(H1) * 2.0 + inside(H2) * 3.0 + inside(H3) * 4.0 + inside(H4) * 5.0;

    const amrex::Real u = 1.0;
    const amrex::Real v = 1.0;
    const amrex::Real w = 1.0;
    const amrex::Real P = 1.0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*rho1, Alpha2*rho2,
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
// one line to <plot directory>/interface.csv with
//   V_2        volume of material 2, int (1 - alpha_1) dV
//   a1_min     min alpha_1
//   a1_max     max alpha_1
//   N_I        interface thickness in finest cells, (1/dx) int H(alpha_1) dV / int |grad alpha_1| dV,
//              with H = 1 where 0.05 <= alpha_1 <= 0.95
// All quantities are taken over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const int ia = INDEX_VolumeFraction1;

    // One ghost cell, for the gradient
    const auto r = Measure(o, {Sum, Sum, Sum, Min, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1 = U(i,j,k,ia);
        AMREX_D_TERM(const Real gx = 0.5 * (U(i+1,j,k,ia) - U(i-1,j,k,ia)) / c.dx[0];,
                     const Real gy = 0.5 * (U(i,j+1,k,ia) - U(i,j-1,k,ia)) / c.dx[1];,
                     const Real gz = 0.5 * (U(i,j,k+1,ia) - U(i,j,k-1,ia)) / c.dx[2];)
        const Real H = (a1 >= 0.05 && a1 <= 0.95) ? 1.0 : 0.0;

        return {(1.0 - a1) * c.dV, H * c.dV, std::sqrt(AMREX_D_TERM(gx*gx, + gy*gy, + gz*gz)) * c.dV,
                a1, a1};
    });

    const Real dx_f = FinestGeom(o).CellSize(0);
    WriteRow(o, "interface.csv", {{"V_2", r[0]}, {"a1_min", r[3]}, {"a1_max", r[4]},
                                  {"N_I", (r[2] > 0.0) ? r[1] / r[2] / dx_f : 0.0}});
}

#endif  // IC_H_
