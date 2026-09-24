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

static_assert(NPHASE == 5, "Advection-Viscous-5Eq-N sets up five materials; build with NPHASE=5.");


// =======================================
// Initial condition
//
// Four shapes in a background of material 1, each filled with its own material: Zalesak's
// slotted disk (Zalesak, JCP 31, 1979) of material 2, a square of material 3, a disk of
// material 4 and an equilateral triangle of material 5. Each shape is a tanh profile of
// width eta = Eta_Multiplier * dx_FinestLevel about its signed-distance function, and the
// volume fractions are normalized to sum to one. The materials have densities 1 to 5 at
// a uniform pressure and move with the uniform velocity (1, 1), so they differ in
// temperature and heat conducts across the interfaces as they are advected.
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

    // Signed distance of each material's shape (positive inside), and its density
    amrex::Array<amrex::Real, NPHASE> Psi;
    const amrex::Array<amrex::Real, NPHASE> rho_phase = {1.0, 2.0, 3.0, 4.0, 5.0};

    // ---------------------------------------
    // Material 2: Zalesak's slotted disk, center (0.4, -0.4), radius 0.4, slot of width
    // 0.125 open at the bottom and reaching 0.25 above the center
    // ---------------------------------------
    amrex::Real X = x - 0.4;
    amrex::Real Y = y + 0.4;
    const amrex::Real r_disk   = 0.4;
    const amrex::Real slot_w   = 0.125;
    const amrex::Real slot_d   = 0.65;
    const amrex::Real psi_disk = r_disk - std::sqrt(X*X + Y*Y);
    const amrex::Real psi_slot = std::min(0.5*slot_w - std::abs(X), (slot_d - r_disk) - Y);
    Psi[1] = std::min(psi_disk, -psi_slot);

    // ---------------------------------------
    // Material 3: square of side 0.4, center (0.5, 0.5)
    // ---------------------------------------
    Psi[2] = 0.2 - std::max(std::abs(x - 0.5), std::abs(y - 0.5));

    // ---------------------------------------
    // Material 4: disk of radius 0.2, center (-0.5, -0.5)
    // ---------------------------------------
    Psi[3] = 0.2 - std::sqrt((x + 0.5)*(x + 0.5) + (y + 0.5)*(y + 0.5) + z*z);

    // ---------------------------------------
    // Material 5: equilateral triangle of side 0.5 pointing in -y, with its upper edge at
    // y = 0.5 + h/3 and its lower vertex at (-0.5, 0.5 - h/3), where h = (sqrt(3)/2) 0.75
    // ---------------------------------------
    X = x + 0.5;
    Y = y - 0.5;
    const amrex::Real h  = std::sqrt(3.0) * 0.5 * 0.75;
    const amrex::Real d1 =  (h/3.0) - Y;
    const amrex::Real d2 =  (std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    const amrex::Real d3 = -(std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    Psi[4] = std::min(d1, std::min(d2, d3));

    // ---------------------------------------
    // Material 1: everything outside the shapes
    // ---------------------------------------
    Psi[0] = -std::max({Psi[1], Psi[2], Psi[3], Psi[4]});

    // ---------------------------------------
    // Volume fractions, normalized to sum to one
    // ---------------------------------------
    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho;
    amrex::Real Sum_Alpha = 0.0;
    for (int n = 0; n < NPHASE; ++n) {
        Alpha[n]   = 0.5*(1.0 + std::tanh(Psi[n]/(std::sqrt(2.0)*eta))*(1.0 - eps));
        Sum_Alpha += Alpha[n];
    }
    for (int n = 0; n < NPHASE; ++n) {
        Alpha[n]    /= Sum_Alpha;
        Alpha_rho[n] = Alpha[n]*rho_phase[n];
    }

    const amrex::Real u = 1.0;
    const amrex::Real v = 1.0;
    const amrex::Real w = 1.0;
    const amrex::Real P = 1.0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/materials.csv with
//   V_1 ... V_N   volume of each material, int alpha_p dV
//   alpha_min     min over the materials of alpha_p
//   sum_err       max |1 - sum_p alpha_p|
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

    const auto V = Measure<NPHASE>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<NPHASE>
    {
        Values<NPHASE> v;
        for (int n = 0; n < NPHASE; ++n) { v[n] = U(c.i,c.j,c.k,ia+n) * c.dV; }
        return v;
    });

    const auto r = Measure(o, {Min, Max}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<2>
    {
        Real lo = 1.0e30, s = 0.0;
        for (int n = 0; n < NPHASE; ++n) {
            lo = std::min(lo, U(c.i,c.j,c.k,ia+n));
            s += U(c.i,c.j,c.k,ia+n);
        }
        return {lo, std::abs(1.0 - s)};
    });

    std::vector<std::pair<std::string, Real>> columns;
    for (int n = 0; n < NPHASE; ++n) { columns.push_back({"V_" + std::to_string(n+1), V[n]}); }
    columns.insert(columns.end(), {{"alpha_min", r[0]}, {"sum_err", r[1]}});
    WriteRow(o, "materials.csv", columns);
}

#endif  // IC_H_
