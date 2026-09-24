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
// Shock-induced instabilities of two successive interfaces (Liang and Luo, JFM 955, 2023;
// section 4.3 of the COMPAS paper). Material 1 (SF6) fills x < x_12 + amp cos(2 pi y/lambda),
// material 2 (an SF6-air mixture) is a layer of thickness L_2, and material 3 (air) is beyond it. A
// shock of Mach number prob.Mach travels in +x through material 1 from x = prob.x_shock.
// The domain is half a wavelength in y, with symmetry at y = 0 and y = lambda/2. With
// NPHASE = 4 the fourth material is initially absent, for the reduction-consistency test.
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

    // ---------------------------------------
    // Interfaces: tanh profiles of width eta about the two perturbed surfaces
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-8;
    const amrex::Real k   = 2.0 * PI / parm.lambda;
    const amrex::Real x12 = parm.x_12 + parm.amp * std::cos(k * y);
    const amrex::Real x23 = parm.x_12 + parm.L_2 + parm.phase_23 * parm.amp * std::cos(k * y);

    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho;
    for (int n = 0; n < NPHASE; ++n) { Alpha[n] = eps; }
    const amrex::Real H1 = 0.5 * (1.0 + std::tanh((x12 - x) / (std::sqrt(2.0) * eta)));
    const amrex::Real H3 = 0.5 * (1.0 + std::tanh((x - x23) / (std::sqrt(2.0) * eta)));
    Alpha[0] = amrex::max(H1, eps);
    Alpha[2] = amrex::max(H3, eps);
    Alpha[1] = amrex::max(1.0 - H1 - H3, eps);
    amrex::Real sum = 0.0;
    for (int n = 0; n < NPHASE; ++n) { sum += Alpha[n]; }
    for (int n = 0; n < NPHASE; ++n) { Alpha[n] /= sum; }

    // ---------------------------------------
    // Gases at rest, and the post-shock state of material 1 behind the shock
    // ---------------------------------------
    amrex::Real P = parm.P0;
    amrex::Real u = 0.0;
    amrex::Array<amrex::Real, NPHASE> rho;
    for (int n = 0; n < NPHASE; ++n) { rho[n] = parm.rho[n]; }
    if (x < parm.x_shock) {
        amrex::Real shock[5];
        StiffenedShockConditions(parm.Mach, 0.0, parm.P0, parm.rho[0], parm.gamma[0], parm.pinf[0], shock);
        for (int n = 0; n < NPHASE; ++n) { rho[n] *= shock[0] / parm.rho[0]; }
        u = shock[1];
        P = shock[2];
    }
    for (int n = 0; n < NPHASE; ++n) { Alpha_rho[n] = Alpha[n] * rho[n]; }
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;

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
// one line to <plot directory>/interfaces.csv with
//   x_12, x_23   mean positions of the two interfaces, where alpha_1 and alpha_3 cross 0.5,
//                averaged between the crest (y = 0) and the trough (y = lambda/2)
//   a_12, a_23   their amplitudes, half the crest-to-trough distance (signed)
//   N_I          thickness of material 2 in cells, (1/dx) int H(alpha_2) dA / int |grad alpha_2| dA,
//                with H = 1 where 0.05 <= alpha_2 <= 0.95
//   sum_err      max |1 - sum_p alpha_p|
//   alpha_last   max alpha of the last material (the absent one when NPHASE = 4)
// Crossings are taken on the finest level; the integrals over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const int  ia  = INDEX_VolumeFraction1;
    const Real big = 1.0e30;

    // One ghost cell, for the neighbors in the crossings and the gradient
    const auto v = Measure(o, {Max, Max, Max, Max, Sum, Sum, Max, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<8>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1 = U(i,j,k,ia), a2 = U(i,j,k,ia+1), a3 = U(i,j,k,ia+2);
        const Real a1e = U(i+1,j,k,ia), a3e = U(i+1,j,k,ia+2);

        // alpha = 0.5 crossings toward the +x neighbor, on the bottom (crest) and top (trough) rows
        Real c12 = -big, c23 = -big;
        if (c.finest) {
            if ((a1 - 0.5) * (a1e - 0.5) <= 0.0 && a1 != a1e) { c12 = c.x[0] + c.dx[0] * (0.5 - a1) / (a1e - a1); }
            if ((a3 - 0.5) * (a3e - 0.5) <= 0.0 && a3 != a3e) { c23 = c.x[0] + c.dx[0] * (0.5 - a3) / (a3e - a3); }
        }
        const bool crest = (j == c.domain.smallEnd(1)), trough = (j == c.domain.bigEnd(1));

        // Thickness of material 2
        const Real gx = 0.5 * (U(i+1,j,k,ia+1) - U(i-1,j,k,ia+1)) / c.dx[0];
        const Real gy = 0.5 * (U(i,j+1,k,ia+1) - U(i,j-1,k,ia+1)) / c.dx[1];
        const Real H  = (a2 >= 0.05 && a2 <= 0.95) ? 1.0 : 0.0;

        Real s = 0.0;
        for (int n = 0; n < NPHASE; ++n) { s += U(i,j,k,ia+n); }

        return {crest ? c12 : -big, trough ? c12 : -big, crest ? c23 : -big, trough ? c23 : -big,
                H * c.dV, std::sqrt(gx*gx + gy*gy) * c.dV, std::abs(1.0 - s), U(i,j,k,ia+NPHASE-1)};
    });

    const Real dx_f = FinestGeom(o).CellSize(0);
    WriteRow(o, "interfaces.csv", {{"x_12", 0.5 * (v[0] + v[1])}, {"x_23", 0.5 * (v[2] + v[3])},
                                   {"a_12", 0.5 * (v[0] - v[1])}, {"a_23", 0.5 * (v[2] - v[3])},
                                   {"N_I", (v[5] > 0.0) ? v[4] / v[5] / dx_f : 0.0},
                                   {"sum_err", v[6]}, {"alpha_last", v[7]}});
}

#endif  // IC_H_
