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
// A pool of heavy fluid 2 fills y < 0 below light fluid 1, and the solid read from the
// STL file (placed in Parm::DynamicInit) is also filled with fluid 2. The volume fraction
// is a tanh profile of width eta about the union of the pool and the solid, built from the
// smoothed signed distance to the mesh. In 2D the plane z = 0 cuts the mesh, and the
// distance inside it is bounded by the distance to its faces in z, so the mesh must be
// several eta thick in z for alpha_2 to reach 1 inside. Everything is at rest, and the
// pressure is hydrostatic for the pool alone, P0 + rho g_y y.
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
    // Volume fraction: fluid 1 where both the distance above the pool surface (y) and
    // the distance outside the solid (-Psi_STL, positive outside) are positive
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;

    const StlGeom::Vector3 P1{float(x), float(y), float(z)};
    const amrex::Real Psi_STL  = amrex::Real(stl.SignedDistanceBVH(P1, eta, true));
    const amrex::Real Psi_Surf = y;
    const amrex::Real Psi      = std::min(-Psi_STL, Psi_Surf);

    const amrex::Real Alpha1 = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // Pool alone, for the hydrostatic pressure
    const amrex::Real Alpha1_Surf = 0.5 * (1.0 + std::tanh(Psi_Surf / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha2_Surf = 1.0 - Alpha1_Surf;

    // ---------------------------------------
    // Material states
    // ---------------------------------------
    const amrex::Real Alpha1_rho1 = Alpha1 * parm.rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * parm.rho2;

    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
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
// User source term: gravity
//
// Adds the body force rho g to the momentum and its work rho g . u to the energy.
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

    amrex::Real rho = Alpha1_rho1 + Alpha2_rho2;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> F_Grav = {AMREX_D_DECL(rho*parm.gx, rho*parm.gy, rho*parm.gz)};
    amrex::Real E_Grav = AMREX_D_TERM(F_Grav[0] * u, + F_Grav[1] * v, + F_Grav[2] * w);

    Residual[0] = 0.0;
    Residual[1] = 0.0;
    AMREX_D_TERM(
    Residual[2] = F_Grav[0];,
    Residual[3] = F_Grav[1];,
    Residual[4] = F_Grav[2];)
    Residual[2 + AMREX_SPACEDIM] = E_Grav;
    Residual[3 + AMREX_SPACEDIM] = 0.0;
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/solid.csv with, for fluid 2 above the initial pool surface
// (y > 0), which is the STL solid until it reaches the pool,
//   V_2    its volume, int alpha_2 dV (an area per unit depth in 2D)
//   y_c    the y-coordinate of its centroid
//   v_c    its mean vertical velocity, int alpha_2 v dV / V_2
// The integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto v = Measure<3>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<3>
    {
        const int  i = c.i, j = c.j, k = c.k;
        if (c.x[1] <= 0.0) { return {0.0, 0.0, 0.0}; }
        const Real a2 = 1.0 - U(i,j,k,INDEX_VolumeFraction1);
        const Real v  = U(i,j,k,INDEX_MomentumY) / (U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2));
        return {a2 * c.dV, a2 * c.x[1] * c.dV, a2 * v * c.dV};
    });
    WriteRow(o, "solid.csv", {{"V_2", v[0]}, {"y_c", (v[0] > 0.0) ? v[1] / v[0] : 0.0},
                              {"v_c", (v[0] > 0.0) ? v[2] / v[0] : 0.0}});
}

#endif  // IC_H_
