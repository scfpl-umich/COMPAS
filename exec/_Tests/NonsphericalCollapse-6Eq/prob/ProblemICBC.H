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
// Near-wall bubble collapse (Johnsen and Colonius, J. Fluid Mech. 629, 2009). A sphere of
// gas (material 1) of radius r0 is centered at (-delta, 0, 0), in units of r0, in water
// (material 2) at rest. The gas is at P_gas, and the liquid pressure rises from P_gas at
// the bubble surface to P_liquid far away as P_liquid + (r0/R)(P_gas - P_liquid), with R
// the distance from the bubble center. With Eta_Multiplier = 0 the interface is sharp.
// Densities and pressures are scaled to the units set in Parm::DynamicInit.
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
    // Volume fraction: tanh profile of width eta about the sphere R = r0
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;
    const amrex::Real R   = std::sqrt((x+parm.delta)*(x+parm.delta) + y*y + z*z);
    const amrex::Real Psi = parm.r0/parm.L0 - R;

    const amrex::Real Alpha1 = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Material states: gas at P_gas, liquid with the 1/R pressure profile, all at rest
    // ---------------------------------------
    const amrex::Real rho1 = parm.rho_gas;
    const amrex::Real rho2 = parm.rho_liquid;

    const amrex::Real P1 = parm.P_gas;
    const amrex::Real P2 = parm.P_liquid + parm.r0/parm.L0/R * (parm.P_gas - parm.P_liquid);

    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;

    // ---------------------------------------
    // Scale, convert to primitive and then conservative
    // ---------------------------------------
    amrex::Real Alpha1_rho1 = Alpha1 * rho1;
    amrex::Real Alpha2_rho2 = Alpha2 * rho2;

    amrex::Real Alpha1_P1 = Alpha1 * P1;
    amrex::Real Alpha2_P2 = Alpha2 * P2;

    Alpha1_rho1 /= parm.rho0; Alpha2_rho2 /= parm.rho0; Alpha1_P1 /= parm.P0; Alpha2_P2 /= parm.P0;
    AMREX_D_TERM(u /= parm.u0;, v /= parm.u0;, w /= parm.u0;)

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1_rho1, Alpha2_rho2,
                    AMREX_D_DECL(u, v, w),
                    Alpha1_P1, Alpha2_P2,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Spatio-temporal refinement tagging
//
// Tags the box x < -delta + 1.25, y < 1.25, z < 1.25 (in units of r0) around the quarter
// bubble on every level, so that the bubble and its surroundings are at the finest level.
// =======================================

#define SPATIO_TEMPORAL_TAG_USER_FUNC_
#define SPACE_TIME_REF
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
spatio_temporal_tag (int & tagged,
                     int lev,
                     amrex::Real x,
                     amrex::Real y,
                     amrex::Real z,
                     amrex::Real t,
                     Parm const& parm)
{
    amrex::ignore_unused(lev, t);
    if (y < 1.25 && z < 1.25 && x < -parm.delta + 1.25){
        tagged = TagBox::SET;
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/bubble.csv with
//   V        gas volume, 4 int alpha_1 dV for the full bubble
//   R_eq     equivalent radius, (3 V / 4 pi)^(1/3)
//   x_c      x-coordinate of the gas centroid, its distance from the wall
//   u_c      gas-averaged x-velocity, int alpha_1 u dV / int alpha_1 dV
//   P_wall   maximum mixture pressure in the cells next to the wall, divided by P_liquid
// in the units of the inputs file (r0, t_c and r0/t_c). The integrals are over the
// uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto v = Measure(o, {Sum, Sum, Sum, Max}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<4>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1 = U(i,j,k,INDEX_VolumeFraction1);
        const Real u  = U(i,j,k,INDEX_MomentumX) / (U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2));

        // Mixture pressure in the cells next to the wall at x = 0
        Real P = -1.0e30;
        if (i == c.domain.smallEnd(0)) {
            amrex::Array<Real, NSTATE> C;
            for (int n = 0; n < NSTATE; ++n) { C[n] = U(i,j,k,n); }
            Conservative2Pressure(P, C, parm);
        }
        return {a1 * c.dV, a1 * c.x[0] * c.dV, a1 * u * c.dV, P};
    });

    const Real V = 4.0 * v[0];
    WriteRow(o, "bubble.csv", {{"V", V}, {"R_eq", std::cbrt(3.0 * V / (4.0 * PI))},
                               {"x_c", (v[0] > 0.0) ? v[1] / v[0] : 0.0},
                               {"u_c", (v[0] > 0.0) ? v[2] / v[0] : 0.0},
                               {"P_wall", v[3] * HostParm(o).P0 / HostParm(o).P_liquid}});
}

#endif  // IC_H_
