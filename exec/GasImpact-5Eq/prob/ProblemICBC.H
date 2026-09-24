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
// and is at twice the ambient pressure. The rest of the domain is material 2 at rest, with
// density prob.rho2 in the plate 1.0 <= x <= 1.1 and prob.rho3 elsewhere. The domain starts
// at y = 1, so only the upper half of the block is simulated. Interfaces are tanh profiles
// of width eta.
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

    const amrex::Real Alpha1 = Alpha_rect;
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Material states: the block moves in +x at twice the ambient pressure
    // ---------------------------------------
    const amrex::Real rho1 = parm.rho1;
    const amrex::Real rho2 = parm.rho2*Alpha_slab + parm.rho3 * (1.0 - Alpha_slab);
    const amrex::Real u = 0.2*Alpha_rect;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    const amrex::Real P = 2.0 * Alpha_rect + (1.0 - Alpha_rect);

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
// one line to <plot directory>/block.csv with
//   V_1   volume of material 1, int alpha_1 dV (area per unit depth in 2D), in the
//         simulated half of the domain
//   x_1   its centroid, int alpha_1 x dV / V_1
//   u_1   its mean x velocity, int alpha_1 rho_1 u dV / int alpha_1 rho_1 dV
// The integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    // V_1, int alpha_1 x dV, the mass of material 1 and int alpha_1 rho_1 u dV
    const auto v = Measure<4>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<4>
    {
        const Real a1   = U(c.i,c.j,c.k,INDEX_VolumeFraction1);
        const Real a1r1 = U(c.i,c.j,c.k,INDEX_Mass1);
        const Real u    = U(c.i,c.j,c.k,INDEX_MomentumX) / (a1r1 + U(c.i,c.j,c.k,INDEX_Mass2));
        return {a1 * c.dV, a1 * c.x[0] * c.dV, a1r1 * c.dV, a1r1 * u * c.dV};
    });
    WriteRow(o, "block.csv", {{"V_1", v[0]}, {"x_1", (v[0] > 0.0) ? v[1] / v[0] : 0.0},
                                             {"u_1", (v[2] > 0.0) ? v[3] / v[2] : 0.0}});
}

#endif  // IC_H_
