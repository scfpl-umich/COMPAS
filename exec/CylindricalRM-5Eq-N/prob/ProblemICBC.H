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
// Three-material cylindrical Richtmyer-Meshkov problem on the quarter plane x, y >= 0.
// A cylindrical shock of Mach number prob.Mach starts at r = prob.Rs in material 1 (air)
// and converges on a heavy shell of material 2 (SF6-like) around a light core of
// material 3 (helium-like). The two interfaces are
//     R_12(theta) = prob.R_12 + prob.a_12 cos(prob.mode_12 theta)   (material 1 / material 2)
//     R_23(theta) = prob.R_23 + prob.a_23 cos(prob.mode_23 theta)   (material 2 / material 3)
// with even mode numbers, so that both axes are lines of symmetry.
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
    // Polar coordinates and interface width
    // ---------------------------------------
    const amrex::Real r     = std::sqrt(x*x + y*y);
    const amrex::Real theta = std::atan2(y, x);
    const amrex::Real eta   = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps   = 1.0e-8;

    // ---------------------------------------
    // Perturbed interfaces: tanh profiles of width eta.
    // H_12 = 1 inside R_12 (materials 2 and 3), H_23 = 1 inside R_23 (material 3)
    // ---------------------------------------
    const amrex::Real R12  = parm.R_12 + parm.a_12 * std::cos(parm.mode_12 * theta);
    const amrex::Real R23  = parm.R_23 + parm.a_23 * std::cos(parm.mode_23 * theta);
    const amrex::Real H_12 = 0.5 * (1.0 + std::tanh((R12 - r) / (std::sqrt(2.0) * eta)));
    const amrex::Real H_23 = 0.5 * (1.0 + std::tanh((R23 - r) / (std::sqrt(2.0) * eta)));

    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho;
    Alpha[0] = amrex::max(1.0 - H_12, eps);
    Alpha[1] = amrex::max(H_12 - H_23, eps);
    Alpha[2] = amrex::max(H_23, eps);
    amrex::Real sum = 0.0;
    for (int n = 0; n < NPHASE; ++n) { sum += Alpha[n]; }
    for (int n = 0; n < NPHASE; ++n) { Alpha[n] /= sum; }

    // ---------------------------------------
    // Materials at rest at P0, and the post-shock state of material 1 outside r = Rs,
    // from the Rankine-Hugoniot relations, moving inward
    // ---------------------------------------
    amrex::Real P   = parm.P0;
    amrex::Real u_r = 0.0;
    amrex::Array<amrex::Real, NPHASE> rho;
    for (int n = 0; n < NPHASE; ++n) { rho[n] = parm.rho[n]; }
    if (r > parm.Rs) {
        amrex::Real shock[5];
        StiffenedShockConditions(parm.Mach, 0.0, parm.P0, parm.rho[0], parm.gamma[0], parm.pinf[0], shock);
        for (int n = 0; n < NPHASE; ++n) { rho[n] *= shock[0] / parm.rho[0]; }
        u_r = -shock[1];
        P   = shock[2];
    }
    for (int n = 0; n < NPHASE; ++n) { Alpha_rho[n] = Alpha[n] * rho[n]; }
    const amrex::Real u = (r > 0.0) ? u_r * x / r : 0.0;
    const amrex::Real v = (r > 0.0) ? u_r * y / r : 0.0;
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
//   R_12, R_23   equivalent radii of the two interfaces, from the areas A = pi R^2 / 4 of
//                materials 2 and 3 together and of material 3 alone
//   h_12, h_23   half the radial extent of the alpha_1 = 0.5 and alpha_3 = 0.5 contours,
//                a measure of the perturbation amplitude and later of the mixing width
// The areas are integrals over the composite grid; the contours are taken on the finest level.
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

    // One ghost cell, for the alpha = 0.5 crossings toward the +x and +y neighbors
    const auto v = Measure(o, {Sum, Sum, Min, Max, Min, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<6>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1 = U(i,j,k,ia), a2 = U(i,j,k,ia+1), a3 = U(i,j,k,ia+2);
        const Real r  = std::sqrt(c.x[0]*c.x[0] + c.x[1]*c.x[1]);

        bool on_12 = false, on_23 = false;
        if (c.finest) {
            const Real s1 = a1 - 0.5, s3 = a3 - 0.5;
            on_12 = (s1 * (U(i+1,j,k,ia)   - 0.5) <= 0.0) || (s1 * (U(i,j+1,k,ia)   - 0.5) <= 0.0);
            on_23 = (s3 * (U(i+1,j,k,ia+2) - 0.5) <= 0.0) || (s3 * (U(i,j+1,k,ia+2) - 0.5) <= 0.0);
        }

        return {(a2 + a3) * c.dV, a3 * c.dV,
                on_12 ? r : big, on_12 ? r : -big,
                on_23 ? r : big, on_23 ? r : -big};
    });

    const Real R_12 = std::sqrt(4.0 * v[0] / PI);
    const Real R_23 = std::sqrt(4.0 * v[1] / PI);
    const Real h_12 = (v[3] > v[2]) ? 0.5 * (v[3] - v[2]) : 0.0;
    const Real h_23 = (v[5] > v[4]) ? 0.5 * (v[5] - v[4]) : 0.0;
    WriteRow(o, "interfaces.csv", {{"R_12", R_12}, {"R_23", R_23}, {"h_12", h_12}, {"h_23", h_23}});
}

#endif  // IC_H_
