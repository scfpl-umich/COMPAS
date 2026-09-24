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
// Collapse of a gas bubble below a free surface. A bubble of unit diameter centered at the
// origin (phase 2) sits in water (phase 1) at the pressure prob.P0, with its gas at P0/30.
// Air (phase 3) at P0 fills y > 1.5, so the free surface is three bubble radii above the
// bubble center. Everything starts at rest. The state is nondimensionalized by prob.rho0,
// prob.P0 and u0 = sqrt(P0/rho0) (see Parm.H).
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
    // Interfaces: tanh profiles of width eta about the bubble surface (radius 0.5) and the
    // free surface (y = 1.5). No floor is applied to the volume fractions.
    // ---------------------------------------
    amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    amrex::Real Psi = 0.5 - std::sqrt(x*x + y*y + z*z);   // > 0 inside the bubble

    amrex::Real Alpha2 = 1.0 - 0.5 * (1.0 + std::tanh(-Psi / (std::sqrt(2.0) * eta)));  // bubble gas
    amrex::Real Alpha3 = 0.5 * (1.0 + std::tanh((y - 1.5) / (std::sqrt(2.0) * eta)));   // air
    amrex::Real Alpha1 = 1.0 - Alpha2 - Alpha3;                                         // water

    // ---------------------------------------
    // Material states (kg/m^3 and Pa), all at rest
    // ---------------------------------------
    amrex::Real rho1 = parm.rho0;
    amrex::Real rho2 = 1.0;
    amrex::Real rho3 = 1.0;

    amrex::Real P1 = parm.P0;
    amrex::Real P2 = parm.P0/30.0;
    amrex::Real P3 = parm.P0;

    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;

    // ---------------------------------------
    // Nondimensionalize, then convert to primitive and conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NPHASE> Alphas     = {Alpha1, Alpha2, Alpha3};
    amrex::Array<amrex::Real, NPHASE> Alpha_rhos = {Alpha1 * rho1 / parm.rho0, Alpha2 * rho2 / parm.rho0, Alpha3 * rho3 / parm.rho0};
    amrex::Array<amrex::Real, NPHASE> Alpha_Ps   = {Alpha1 * P1 / parm.P0, Alpha2 * P2 / parm.P0, Alpha3 * P3 / parm.P0};
    AMREX_D_TERM(u /= parm.u0;, v /= parm.u0;, w /= parm.u0;)

    amrex::Real P = Alpha_Ps[0] + Alpha_Ps[1] + Alpha_Ps[2];

    amrex::Array<amrex::Real, NSTATE> Primitive;
#if (PHYSICS==SIXEQS_IE_NPHASE)
    State2Primitive(Primitive,
                    Alphas,
                    Alpha_rhos,
                    Alpha_Ps,
       AMREX_D_DECL(u, v, w),
                    P);
#elif (PHYSICS==FIVEEQS_NPHASE)
    State2Primitive(Primitive,
                    Alpha_rhos,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alphas);
#endif
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Optional: spatio-temporal refinement tagging
//
// Every level is refined within a distance 2 of the bubble center, which covers the bubble
// and the part of the free surface above it.
// =======================================

#define SPATIO_TEMPORAL_TAG_USER_FUNC_
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void spatio_temporal_tag(int& TagVal,
                         int lev,
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z,
                         amrex::Real t,
                         Parm const& parm)
{
    amrex::ignore_unused(lev, t, parm);

    amrex::Real r = std::sqrt(x*x + y*y + z*z);
    if (r < 2.0){
        TagVal = TagBox::SET;
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/bubble.csv with, in the nondimensional units of the case,
//   V_b      amount of bubble gas, int alpha_2 dV (the area of the half bubble in 2D)
//   y_c      y centroid of the bubble gas, int alpha_2 y dV / V_b
//   y_fs     height of the free surface on the symmetry axis, where alpha_3 crosses 0.5
//   P_max    largest mixture pressure in the domain
// The integrals and extrema are taken over the composite grid.
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

    // One ghost cell, for the +y neighbor in the free-surface crossing
    const auto v = Measure(o, {Sum, Sum, Max, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<4>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a2 = U(i,j,k,ia+1), a3 = U(i,j,k,ia+2), a3n = U(i,j+1,k,ia+2);

        // Free surface: alpha_3 = 0.5 crossing toward the +y neighbor, in the cells next to the symmetry axis
        bool axis = (i == c.domain.smallEnd(0));
#if (AMREX_SPACEDIM == 3)
        axis = axis && (k == c.domain.smallEnd(2));
#endif
        Real y_fs = -big;
        if (axis && (a3 - 0.5) * (a3n - 0.5) <= 0.0 && a3 != a3n) {
            y_fs = c.x[1] + c.dx[1] * (0.5 - a3) / (a3n - a3);
        }

        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(i,j,k,n); }
        Real P;
        Conservative2Pressure(P, Cons, parm);

        return {a2 * c.dV, a2 * c.x[1] * c.dV, y_fs, P};
    });

    WriteRow(o, "bubble.csv", {{"V_b", v[0]}, {"y_c", (v[0] > 0.0) ? v[1] / v[0] : 0.0},
                               {"y_fs", (v[2] > -big) ? v[2] : 0.0}, {"P_max", v[3]}});
}

#endif  // IC_H_
