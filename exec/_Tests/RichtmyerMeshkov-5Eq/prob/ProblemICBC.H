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
// Single-mode Richtmyer-Meshkov instability. Material 1 (light) fills y > a0 cos(2 pi x/lambda),
// with lambda = 1, and material 2 (heavy) is below it. Both are at rest at unit pressure.
// A shock of Mach number prob.Mach starts at y = 1.1 a0 and travels in -y through material 1,
// so the post-shock state of material 1 is set above it. The interface is a tanh profile of
// width eta = PhaseField.Eta_Multiplier finest cells.
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
    // Interface: tanh profile about y = a0 cos(2 pi x/lambda)
    // ---------------------------------------
    const amrex::Real eta    = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps    = 1.0e-12;
    const amrex::Real lambda = 1.0;
    const amrex::Real Psi    = y - parm.a0 * std::cos(2.0 * PI * x / lambda);

    const amrex::Real Alpha1 = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Gases at rest, and the post-shock state of material 1 above the shock
    // ---------------------------------------
    const amrex::Real P_Unshocked = 1.0;

    amrex::Real rho1 = parm.rho1;
    amrex::Real rho2 = parm.rho2;
    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;
    amrex::Real P = P_Unshocked;

    const amrex::Real y_Shock = 1.1 * parm.a0;
    if (y > y_Shock) {
        amrex::Real shock[5];
        StiffenedShockConditions(parm.Mach, 0.0, P_Unshocked, parm.rho1, parm.gamma_1, parm.pinf_1, shock);
        rho1 = shock[0];
        v    = -shock[1];
        P    = shock[2];
    }

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    const amrex::Real Alpha1_rho1 = Alpha1 * rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * rho2;

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1_rho1, Alpha2_rho2,
                    AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Pressure of cell (i,j,k), for the diagnostics below
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real CellPressure(amrex::Array4<const amrex::Real> const& U, int i, int j, int k, Parm const& parm)
{
    amrex::Array<amrex::Real, NSTATE> Conservative, Primitive;
    for (int n = 0; n < NSTATE; ++n) { Conservative[n] = U(i,j,k,n); }
    Conservative2Primitive(Primitive, Conservative, parm);
    return Primitive[INDEX_Energy];
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps, and otherwise with
// each plotfile. Each call appends one line to <plot directory>/out1D_<case name>.csv with
//   I_L, I_R        interface position, where alpha_1 crosses 0.5 in y, at the crest (x = 0)
//                   and at the trough (x = lambda/2), taken on the finest level
//   A               amplitude, |I_L - I_R| / 2
//   Gamma, Gamma_F  circulation, int omega dA with omega = dv/dx - du/dy, over the domain and
//                   over the interface region (alpha_1 alpha_2 > 1e-5)
//   Ens, Ens_F      enstrophy, int omega^2/2 dA, over the same two regions
//   L               interface length, int |grad alpha_1| dA
//   TV_P, TV_P_F    total variation of pressure, int |grad P| dA, over the same two regions
// The integrals are over the composite grid.
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

    // One ghost cell, for the neighbors in the crossings and the gradients
    const auto v = Measure(o, {Max, Max, Sum, Sum, Sum, Sum, Sum, Sum, Sum}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<9>
    {
        const int i = c.i, j = c.j, k = c.k;

        // alpha_1 = 0.5 crossing toward the -y neighbor, in the crest and trough columns
        Real cross = -big;
        const Real aL = U(i,j-1,k,ia), aR = U(i,j,k,ia);
        const int  klo = (AMREX_SPACEDIM == 3) ? c.domain.smallEnd(AMREX_SPACEDIM-1) : 0;
        if (c.finest && k == klo && aL <= 0.5 && aR >= 0.5 && aR > aL) {
            cross = c.x[1] - c.dx[1] + c.dx[1] * (0.5 - aL) / (aR - aL);
        }

        // Vorticity
        const Real rho_xm = U(i-1,j,k,INDEX_Mass1) + U(i-1,j,k,INDEX_Mass2);
        const Real rho_xp = U(i+1,j,k,INDEX_Mass1) + U(i+1,j,k,INDEX_Mass2);
        const Real rho_ym = U(i,j-1,k,INDEX_Mass1) + U(i,j-1,k,INDEX_Mass2);
        const Real rho_yp = U(i,j+1,k,INDEX_Mass1) + U(i,j+1,k,INDEX_Mass2);
        const Real dvdx = 0.5 * (U(i+1,j,k,INDEX_MomentumY) / rho_xp - U(i-1,j,k,INDEX_MomentumY) / rho_xm) / c.dx[0];
        const Real dudy = 0.5 * (U(i,j+1,k,INDEX_MomentumX) / rho_yp - U(i,j-1,k,INDEX_MomentumX) / rho_ym) / c.dx[1];
        const Real omega = dvdx - dudy;

        // Gradients of alpha_1 and of the pressure
        Real ga2 = 0.0, gp2 = 0.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int di = (d == 0), dj = (d == 1), dk = (d == 2);
            const Real ga = 0.5 * (U(i+di,j+dj,k+dk,ia) - U(i-di,j-dj,k-dk,ia)) / c.dx[d];
            const Real gp = 0.5 * (CellPressure(U, i+di, j+dj, k+dk, parm) - CellPressure(U, i-di, j-dj, k-dk, parm)) / c.dx[d];
            ga2 += ga * ga;
            gp2 += gp * gp;
        }

        const Real a1  = U(i,j,k,ia);
        const Real F   = (a1 * (1.0 - a1) > 1.0e-5) ? 1.0 : 0.0;
        const Real ens = 0.5 * omega * omega;
        const Real tvp = std::sqrt(gp2);

        return {(i == c.domain.smallEnd(0)) ? cross : -big, (i == c.domain.bigEnd(0)) ? cross : -big,
                omega * c.dV, F * omega * c.dV, ens * c.dV, F * ens * c.dV,
                std::sqrt(ga2) * c.dV, tvp * c.dV, F * tvp * c.dV};
    });

    WriteRow(o, "out1D_" + o.case_name + ".csv",
             {{"I_L", v[0]}, {"I_R", v[1]}, {"A", 0.5 * std::abs(v[0] - v[1])},
              {"Gamma", v[2]}, {"Gamma_F", v[3]}, {"Ens", v[4]}, {"Ens_F", v[5]},
              {"L", v[6]}, {"TV_P", v[7]}, {"TV_P_F", v[8]}});
}

#endif  // IC_H_
