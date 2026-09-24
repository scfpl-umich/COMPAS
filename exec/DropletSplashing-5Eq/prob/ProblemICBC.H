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
// A pool of material 2 fills y < a0 cos(2 pi x) under material 1, and a circular droplet of
// material 2 (center (xc, yc), radius R) is at rest above it. The pressure is hydrostatic
// about the pool surface, P = P0 + rho g_y y with the density of the pool or of the gas
// above it, so the droplet starts in the hydrostatic field of the surrounding gas.
// Interfaces are tanh profiles of width eta.
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
    // Pool surface, of wavelength 1, and its volume fraction alone (for the pressure)
    // ---------------------------------------
    const amrex::Real lambda = 1.0;
    amrex::Real Psi = y - parm.a0 * std::cos(2.0 * PI * x / lambda);

    const amrex::Real Alpha1_Surf = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2_Surf = 1.0 - Alpha1_Surf;

    // ---------------------------------------
    // Add the droplet: union of the pool and the circle (level sets negative in material 2)
    // ---------------------------------------
    const amrex::Real Psi_circle = std::sqrt((x - parm.xc)*(x - parm.xc) + (y - parm.yc)*(y - parm.yc)) - parm.R;
    Psi = std::min(Psi, Psi_circle);

    const amrex::Real Alpha1 = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Material states: at rest, hydrostatic pressure
    // ---------------------------------------
    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;

    const amrex::Real Alpha1_rho1 = Alpha1 * parm.rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * parm.rho2;

    const amrex::Real P = parm.P0 + (Alpha1_Surf * parm.rho1 + Alpha2_Surf * parm.rho2) * parm.gy * y;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
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
// Source term: gravity
//
// Adds the body force rho g to the momentum and its work rho g . u to the energy, with
// g = (prob.gx, prob.gy, prob.gz).
// =======================================

#define USER_SOURCE_TERM
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void user_source_term (
    amrex::Array<amrex::Real, NSTATE>& Residual,
    amrex::Array4<const amrex::Real> const& statein,
    int i,
    int j,
    int k,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z,
    amrex::Real t,
    amrex::Real dt,
    AMREX_D_DECL(amrex::Real dx,
                 amrex::Real dy,
                 amrex::Real dz),
    Parm const& parm)
{
    amrex::ignore_unused(x, y, z, t, dt, AMREX_D_DECL(dx, dy, dz));

    amrex::Array<amrex::Real, NSTATE> Conservative;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = statein(i,j,k,iState);
    }
    Conservative2Primitive(Primitive,Conservative,parm);

    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;)
    amrex::Real p;
    amrex::Real Alpha1;
    Primitive2State(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u,v,w), p, Alpha1, Primitive);

    const amrex::Real rho = Alpha1_rho1 + Alpha2_rho2;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> F_Grav = {AMREX_D_DECL(rho*parm.gx, rho*parm.gy, rho*parm.gz)};
    const amrex::Real E_Grav = AMREX_D_TERM(F_Grav[0] * u, + F_Grav[1] * v, + F_Grav[2] * w);

    Residual[INDEX_Mass1] = 0.0;
    Residual[INDEX_Mass2] = 0.0;
    AMREX_D_TERM(
    Residual[INDEX_MomentumX] = F_Grav[0];,
    Residual[INDEX_MomentumY] = F_Grav[1];,
    Residual[INDEX_MomentumZ] = F_Grav[2];)
    Residual[INDEX_Energy]          = E_Grav;
    Residual[INDEX_VolumeFraction1] = 0.0;
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/droplet.csv with integrals of material 2 over the region
// y > |a0| + R/2, which lies clear of the diffuse pool surface and holds the droplet
// and, after impact, the splash:
//   V_2   volume, int alpha_2 dV (area per unit depth in 2D)
//   y_2   centroid height, int alpha_2 y dV / V_2
//   v_2   mean vertical velocity, int alpha_2 rho_2 v dV / int alpha_2 rho_2 dV
// The integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real y_cut = std::abs(HostParm(o).a0) + 0.5 * HostParm(o).R;

    // Above y_cut: V_2, int alpha_2 y dV, the mass of material 2 and int alpha_2 rho_2 v dV
    const auto v = Measure<4>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<4>
    {
        const Real y = c.x[1];
        if (y <= y_cut) { return {0.0, 0.0, 0.0, 0.0}; }
        const Real a2   = 1.0 - U(c.i,c.j,c.k,INDEX_VolumeFraction1);
        const Real a2r2 = U(c.i,c.j,c.k,INDEX_Mass2);
        const Real vy   = U(c.i,c.j,c.k,INDEX_MomentumY) / (U(c.i,c.j,c.k,INDEX_Mass1) + a2r2);
        return {a2 * c.dV, a2 * y * c.dV, a2r2 * c.dV, a2r2 * vy * c.dV};
    });
    WriteRow(o, "droplet.csv", {{"V_2", v[0]}, {"y_2", (v[0] > 0.0) ? v[1] / v[0] : 0.0},
                                               {"v_2", (v[2] > 0.0) ? v[3] / v[2] : 0.0}});
}

#endif  // IC_H_
