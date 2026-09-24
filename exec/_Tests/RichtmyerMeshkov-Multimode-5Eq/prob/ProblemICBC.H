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
// Multimode Richtmyer-Meshkov instability. Material 1 (light) fills
// y > sum_m a_m cos(k_m x + phi_m), and material 2 (heavy) is below it. Both are at rest at
// unit pressure. A shock of Mach number prob.Mach starts at y = 1.1 a0 and travels in -y
// through material 1, so the post-shock state of material 1 is set above it. The modes are
// read at run time in Parm::DynamicInit and read here through device pointers. The interface
// is a tanh profile of width eta = PhaseField.Eta_Multiplier finest cells.
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
    // Interface: tanh profile about the sum of the modes
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;

    amrex::Real interface = 0.0;
    for (int m = 0; m < parm.n_modes; ++m) {
        interface += parm.a[m] * std::cos(parm.kx[m] * x + parm.phase[m]);
    }
    const amrex::Real Psi = y - interface;

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
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps, and otherwise with
// each plotfile. Each call appends one line to <plot directory>/out1D_<case name>.csv with
//   y_min, y_max    lowest and highest cell centers of the mixing zone, 0.01 < alpha_1 < 0.99
//   W               mixing-zone width, y_max - y_min
//   Cy              y centroid of the mixing zone
//   Gamma, Gamma_F  circulation, int omega dA with omega = dv/dx - du/dy, over the domain and
//                   over the interface region (alpha_1 alpha_2 > 1e-5)
//   Ens, Ens_F      enstrophy, int omega^2/2 dA, over the same two regions
// All quantities are over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // One ghost cell, for the neighbors in the vorticity
    const auto v = Measure(o, {Min, Max, Sum, Sum, Sum, Sum, Sum, Sum}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<8>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real a1  = U(i,j,k,INDEX_VolumeFraction1);
        const Real y   = c.x[1];
        const bool mix = (a1 > 0.01 && a1 < 0.99);
        const Real F   = (a1 * (1.0 - a1) > 1.0e-5) ? 1.0 : 0.0;

        // Vorticity
        const Real rho_xm = U(i-1,j,k,INDEX_Mass1) + U(i-1,j,k,INDEX_Mass2);
        const Real rho_xp = U(i+1,j,k,INDEX_Mass1) + U(i+1,j,k,INDEX_Mass2);
        const Real rho_ym = U(i,j-1,k,INDEX_Mass1) + U(i,j-1,k,INDEX_Mass2);
        const Real rho_yp = U(i,j+1,k,INDEX_Mass1) + U(i,j+1,k,INDEX_Mass2);
        const Real dvdx = 0.5 * (U(i+1,j,k,INDEX_MomentumY) / rho_xp - U(i-1,j,k,INDEX_MomentumY) / rho_xm) / c.dx[0];
        const Real dudy = 0.5 * (U(i,j+1,k,INDEX_MomentumX) / rho_yp - U(i,j-1,k,INDEX_MomentumX) / rho_ym) / c.dx[1];
        const Real omega = dvdx - dudy;
        const Real ens   = 0.5 * omega * omega;

        return {mix ? y : big, mix ? y : -big, mix ? y * c.dV : 0.0, mix ? c.dV : 0.0,
                omega * c.dV, F * omega * c.dV, ens * c.dV, F * ens * c.dV};
    });

    WriteRow(o, "out1D_" + o.case_name + ".csv",
             {{"y_min", v[0]}, {"y_max", v[1]}, {"W", v[1] - v[0]}, {"Cy", (v[3] > 0.0) ? v[2] / v[3] : 0.0},
              {"Gamma", v[4]}, {"Gamma_F", v[5]}, {"Ens", v[6]}, {"Ens_F", v[7]}});
}

#endif  // IC_H_
