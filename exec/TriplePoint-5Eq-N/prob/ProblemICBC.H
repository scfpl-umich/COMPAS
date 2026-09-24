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


// Signed distance to the rectangle [xmin, xmax] x [ymin, ymax], negative inside
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
boxLevelSet2D( amrex::Real x, 
               amrex::Real y, 
               amrex::Real xmin,
               amrex::Real ymin,
               amrex::Real xmax, 
               amrex::Real ymax) 
{
    amrex::Real cx = 0.5 * (xmin + xmax);
    amrex::Real cy = 0.5 * (ymin + ymax);
    amrex::Real hx = 0.5 * (xmax - xmin);
    amrex::Real hy = 0.5 * (ymax - ymin);

    amrex::Real dx = std::abs(x - cx) - hx;
    amrex::Real dy = std::abs(y - cy) - hy;

    amrex::Real outside = std::sqrt(std::max(dx, 0.0) * std::max(dx, 0.0) +
                               std::max(dy, 0.0) * std::max(dy, 0.0));

    amrex::Real inside = std::min(std::max(dx, dy), 0.0);

    return outside + inside;
}

// =======================================
// Initial condition: triple-point problem (Galera, Maire and Breil, JCP 229, 2010)
//
// Three gases at rest in [0,7] x [0,3]: a high-pressure driver, material 2 (rho = 1,
// p = 1), fills x < 1; material 1 (rho = 1, p = 0.1) fills x > 1, y < 1.5; material 3
// (rho = 0.125, p = 0.1) fills x > 1, y > 1.5. The different impedances of materials 1
// and 3 shear the driver-gas front and roll up a vortex about the point where the three
// materials meet. The interfaces lie on cell faces and PhaseField.Eta_Multiplier = 0 makes
// them sharp; the quadrature projection of the initial condition gives exact cell averages.
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

    // ---------------------------------------
    // Volume fractions: tanh profiles of the rectangle level sets (positive inside),
    // with material 3 taking the remainder
    // ---------------------------------------
    amrex::Array<amrex::Real, NPHASE> Psi_phase   = {0.0, 0.0, 0.0};
    amrex::Array<amrex::Real, NPHASE> Alpha_phase = {0.0, 0.0, 0.0};
    amrex::Array<amrex::Real, NPHASE> rho_phase   = {1.0, 1.0, 0.125};
    amrex::Array<amrex::Real, NPHASE> p0          = {0.1, 1.0, 0.1};

    Psi_phase[0] = -boxLevelSet2D(x, y,  1.0, -1.0, 8.0, 1.5);
    Psi_phase[1] = -boxLevelSet2D(x, y, -1.0, -1.0, 1.0, 4.0);

    Alpha_phase[0] = 0.5*(1.0 + std::tanh( Psi_phase[0]/(std::sqrt(2.0)*eta) )*(1.0 - eps) );
    Alpha_phase[1] = 0.5*(1.0 + std::tanh( Psi_phase[1]/(std::sqrt(2.0)*eta) )*(1.0 - eps) );
    Alpha_phase[2] = 1 - Alpha_phase[0] - Alpha_phase[1];

    amrex::Real Sum_Alpha = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Sum_Alpha += Alpha_phase[i_Phase];
    }
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_phase[i_Phase] /= Sum_Alpha;
    }

    // ---------------------------------------
    // Material states: all at rest, pressure mixed by volume fraction
    // ---------------------------------------
    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;

    amrex::Real P = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        P += Alpha_phase[i_Phase] * p0[i_Phase];
    }

    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha_phase[i_Phase]*rho_phase[i_Phase];
    }

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha_phase);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/triple_point.csv with
//   V_1, V_2, V_3   volume of each material, int alpha_p dV (area per unit depth in 2D)
//   x_bot, x_top    extent of the driver gas (material 2) along the bottom (y = 0) and
//                   top (y = 3) walls, int alpha_2 dx along the row of cells at each wall
// The integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    // Cross-section normal to x and y (1 in 2D), to turn a wall-row integral into a length
    const Real L_z = (AMREX_SPACEDIM == 3) ? o.Geom(0).ProbLength(AMREX_SPACEDIM-1) : 1.0;

    const auto v = Measure<5>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<5>
    {
        const Real a1 = U(c.i,c.j,c.k,INDEX_VolumeFraction1);
        const Real a2 = U(c.i,c.j,c.k,INDEX_VolumeFraction1+1);
        const Real a3 = U(c.i,c.j,c.k,INDEX_VolumeFraction1+2);
        const Real dL = c.dV / c.dx[1] / L_z;
        return {a1 * c.dV, a2 * c.dV, a3 * c.dV,
                (c.j == c.domain.smallEnd(1)) ? a2 * dL : 0.0,
                (c.j == c.domain.bigEnd(1))   ? a2 * dL : 0.0};
    });
    WriteRow(o, "triple_point.csv", {{"V_1", v[0]}, {"V_2", v[1]}, {"V_3", v[2]}, {"x_bot", v[3]}, {"x_top", v[4]}});
}

#endif  // IC_H_
