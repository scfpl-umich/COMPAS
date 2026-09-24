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
// Material 1 fills the block 0.25 <= x <= 0.75, 0.9 <= y <= 1.1, which moves at u = 0.2
// and is at twice the ambient pressure. Material 2 fills the plate 1.0 <= x <= 1.1, and
// material 3 is the surrounding gas, both at rest. The domain starts at y = 1, so only the
// upper half of the block is simulated. Interfaces are tanh profiles of width eta.
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

    // ---------------------------------------
    // Level sets, positive inside: the block (distance to the nearest edge) and the plate
    // ---------------------------------------
    const amrex::Real xmin = 0.25;
    const amrex::Real xmax = 0.75;
    const amrex::Real ymin = 0.9;
    const amrex::Real ymax = 1.1;
    const amrex::Real Psi_rect = std::min({x - xmin, xmax - x, y - ymin, ymax - y});
    const amrex::Real Psi_slab = 0.5*(0.1) - std::abs(x - (1.0) - 0.5*(0.1));

    const amrex::Real Alpha_rect = 0.5 * (1.0 + std::tanh(Psi_rect / (std::sqrt(2.0)*eta)));
    const amrex::Real Alpha_slab = 0.5 * (1.0 + std::tanh(Psi_slab / (std::sqrt(2.0)*eta)));

    amrex::Array<amrex::Real, NPHASE> Alpha_phase;
    Alpha_phase[0] = Alpha_rect;
    Alpha_phase[1] = Alpha_slab;
    Alpha_phase[2] = 1.0 - Alpha_phase[0] - Alpha_phase[1];

    amrex::Real Sum_Alpha = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Sum_Alpha += Alpha_phase[i_Phase];
    }
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_phase[i_Phase] /= Sum_Alpha;
    }

    // ---------------------------------------
    // Material states: the block moves in +x at twice the ambient pressure
    // ---------------------------------------
    const amrex::Real u = 0.2*Alpha_rect;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    const amrex::Real P = 1.0 + 1.0 * Alpha_phase[0];

    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha_phase[i_Phase]*parm.rho[i_Phase];
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
// one line to <plot directory>/block_plate.csv with, for the block (p = 1) and the plate
// (p = 2),
//   V_p   volume of material p, int alpha_p dV (area per unit depth in 2D), in the
//         simulated half of the domain
//   x_p   its centroid, int alpha_p x dV / V_p
//   u_p   its mean x velocity, int alpha_p rho_p u dV / int alpha_p rho_p dV
// The integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    // For p = 1, 2: V_p, int alpha_p x dV, the mass of material p and int alpha_p rho_p u dV
    const auto v = Measure<8>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<8>
    {
        Real rho = 0.0;
        for (int n = 0; n < NPHASE; ++n) { rho += U(c.i,c.j,c.k,INDEX_Mass1+n); }
        const Real u = U(c.i,c.j,c.k,INDEX_MomentumX) / rho;
        Values<8> s;
        for (int p = 0; p < 2; ++p) {
            const Real a = U(c.i,c.j,c.k,INDEX_VolumeFraction1+p);
            const Real m = U(c.i,c.j,c.k,INDEX_Mass1+p);
            s[4*p] = a * c.dV;  s[4*p+1] = a * c.x[0] * c.dV;  s[4*p+2] = m * c.dV;  s[4*p+3] = m * u * c.dV;
        }
        return s;
    });
    auto ratio = [] (Real a, Real b) { return (b > 0.0) ? a / b : 0.0; };
    WriteRow(o, "block_plate.csv", {{"V_1", v[0]}, {"x_1", ratio(v[1], v[0])}, {"u_1", ratio(v[3], v[2])},
                                    {"V_2", v[4]}, {"x_2", ratio(v[5], v[4])}, {"u_2", ratio(v[7], v[6])}});
}

#endif  // IC_H_
