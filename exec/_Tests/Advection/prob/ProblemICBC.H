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


// One ProblemICBC.H for the four builds of this case (BUILD in the GNUmakefile). The two-phase
// models (5Eq, 6Eq) advect four shapes of material 2 in material 1; the N-phase models (5Eq-N,
// 6Eq-N, NPHASE = 5) advect the four shapes as four materials in material 1. Advection-THINC-5Eq
// and Advection-PF-5Eq run the two-phase problem on the 5Eq build with their own inputs files.

// Uniform pressure and velocity of the flow, shared by the initial condition and the diagnostics
namespace Uniform {
    constexpr amrex::Real P = 1.0;
    constexpr amrex::Real u = 1.0;
    constexpr amrex::Real v = 1.0;
}


#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)

// #######################################
// Two-phase models (5Eq, 6Eq): four shapes of material 2
// #######################################

// =======================================
// Initial condition
//
// Four shapes of material 2 in material 1, each a tanh profile of width eta about a signed
// distance psi (positive inside): a slotted disk (Zalesak, JCP 31, 1979), a square, a circle
// and an equilateral triangle, apex down. Material 1 has density 1 and material 2 has density
// 2, 3, 4 and 5 in the four shapes. The pressure and velocity are uniform, so the interfaces
// only translate (Abgrall, JCP 125, 1996). The volume fractions are kept eps/2 away from 0
// and 1.
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

    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;
    auto profile = [&](amrex::Real psi)
    {
        return 0.5*(1.0 + std::tanh( psi/(std::sqrt(2.0)*eta) )*(1.0 - eps) );
    };

    // Signed distance, volume fraction and density of material 2 in each shape
    amrex::Array<amrex::Real, 4> Psi, Alpha_shape, rho_shape;
    amrex::Real xc, yc, zc, r, L, X, Y;

    // ---------------------------------------
    // Slotted disk: radius 0.4, slot of width 0.125 cut up from the bottom to 0.65 deep
    // ---------------------------------------
    xc = 0.4;
    yc = -0.4;
    r  = 0.4;
    const amrex::Real slot_w = 0.125;
    const amrex::Real slot_d = 0.65;
    X = x - xc;
    Y = y - yc;
    const amrex::Real psi_disk   = r - std::sqrt(X*X + Y*Y);
    const amrex::Real psi_slot_x = (0.5 * slot_w) - std::abs(X);
    const amrex::Real psi_slot_y = (-r + slot_d) - Y;
    const amrex::Real psi_slot   = std::min(psi_slot_x, psi_slot_y);   // positive inside the slot
    Psi[0]       = std::min(psi_disk, -psi_slot);
    rho_shape[0] = 2.0;

    // ---------------------------------------
    // Square of side 0.4
    // ---------------------------------------
    xc = 0.4;
    yc = 0.5;
    L  = 0.4;
    Psi[1]       = 0.5 * L - std::max(std::abs(x - xc), std::abs(y - yc));
    rho_shape[1] = 3.0;

    // ---------------------------------------
    // Circle of radius 0.2 (a sphere in 3D)
    // ---------------------------------------
    xc = -0.5;
    yc = -0.5;
    zc =  0.0;
    r  =  0.2;
    Psi[2]       = r - std::sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
    rho_shape[2] = 4.0;

    // ---------------------------------------
    // Equilateral triangle, apex down: d1, d2 and d3 are the signed distances to its three
    // edges, which gives a height 2h/3 and a side 2L/3
    // ---------------------------------------
    xc = -0.5;
    yc =  0.5;
    L  = 0.75;
    const amrex::Real h = std::sqrt(3.0) * 0.5 * L;
    X = x - xc;
    Y = y - yc;
    const amrex::Real d1 =  (h/3.0) - Y;
    const amrex::Real d2 =  (std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    const amrex::Real d3 = -(std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    Psi[3]       = std::min(d1, std::min(d2, d3));
    rho_shape[3] = 5.0;

    // ---------------------------------------
    // Mixture: material 2 takes the density of the shape it is near (alpha > 1e-5)
    // ---------------------------------------
    for (int n = 0; n < 4; ++n) { Alpha_shape[n] = profile(Psi[n]); }
    auto near = [](amrex::Real alpha) { return (alpha > 1.0e-5) ? 1.0 : 0.0; };

    // The floors of the other three shapes would push Alpha2 past 1 inside a shape
    const amrex::Real Alpha2 = std::min(Alpha_shape[0] + Alpha_shape[1] + Alpha_shape[2] + Alpha_shape[3],
                                        1.0 - 0.5*eps);
    const amrex::Real Alpha1 = 1.0 - Alpha2;
    const amrex::Real rho1   = 1.0;
    const amrex::Real rho2   = near(Alpha_shape[0]) * rho_shape[0] +
                               near(Alpha_shape[1]) * rho_shape[1] +
                               near(Alpha_shape[2]) * rho_shape[2] +
                               near(Alpha_shape[3]) * rho_shape[3];

    const amrex::Real Alpha1_rho1 = Alpha1 * rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * rho2;

#if (PHYSICS == FIVEEQS)
    const amrex::Real u = Uniform::u;
    const amrex::Real v = Uniform::v;
    const amrex::Real w = 0.0;
    const amrex::Real P = Uniform::P;
#else
    const amrex::Real u  = Uniform::u;
    const amrex::Real v  = Uniform::v;
    const amrex::Real w  = 0.0;
    const amrex::Real P1 = Uniform::P;
    const amrex::Real P2 = Uniform::P;
#endif

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
#if (PHYSICS == SIXEQS)
    const amrex::Real Alpha1_P1 = Alpha1 * P1;
    const amrex::Real Alpha2_P2 = Alpha2 * P2;

#endif
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1_rho1, Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
#if (PHYSICS == FIVEEQS)
                    P,
#else
                    Alpha1_P1, Alpha2_P2,
#endif
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/advection.csv with
//   V_2                   volume (area in 2D) of material 2, int (1 - alpha_1) dV
//   alpha_min, alpha_max  extrema of alpha_1
//   dP                    max |P - P0|, with the mixture pressure P = alpha_1 P_1 + alpha_2 P_2
//                         in SIXEQS
//   du                    max |u - u0| and |v - v0|
// where P0, u0 and v0 are the uniform state of the initial condition, which exact advection
// keeps, and, with prob.interface_thickness = 1 (Advection-5Eq and the THINC and Phase-Field
// demos),
//   N_I                   interface thickness in finest cells,
//                         (1/dx) int H(alpha_1) dV / int |grad alpha_1| dV, with H = 1 where
//                         0.05 <= alpha_1 <= 0.95 and the gradient by central differences
// The integrals and the extrema are taken over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto r = Measure(o, {Sum, Min, Max, Max, Max}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Conservative2Primitive(Prim, Cons, parm);
#if (PHYSICS == FIVEEQS)
        Real a1r1, a2r2, P, a1;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(a1r1, a2r2, AMREX_D_DECL(u, v, w), P, a1, Prim);
#else
        Real a1r1, a2r2, a1P1, a2P2, a1;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(a1r1, a2r2, AMREX_D_DECL(u, v, w), a1P1, a2P2, a1, Prim);
        const Real P = a1P1 + a2P2;
#endif

        return {(1.0 - a1) * c.dV, a1, a1, std::abs(P - Uniform::P),
                std::max(std::abs(u - Uniform::u), std::abs(v - Uniform::v))};
    });

    std::vector<std::pair<std::string, Real>> columns{
        {"V_2", r[0]}, {"alpha_min", r[1]}, {"alpha_max", r[2]}, {"dP", r[3]}, {"du", r[4]}};

    if (HostParm(o).interface_thickness) {
        const int ia = INDEX_VolumeFraction1;

        // One ghost cell, for the gradient
        const auto t = Measure(o, {Sum, Sum}, 1,
        [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<2>
        {
            const int  i = c.i, j = c.j, k = c.k;
            const Real a1 = U(i,j,k,ia);
            AMREX_D_TERM(const Real gx = 0.5 * (U(i+1,j,k,ia) - U(i-1,j,k,ia)) / c.dx[0];,
                         const Real gy = 0.5 * (U(i,j+1,k,ia) - U(i,j-1,k,ia)) / c.dx[1];,
                         const Real gz = 0.5 * (U(i,j,k+1,ia) - U(i,j,k-1,ia)) / c.dx[2];)
            const Real H = (a1 >= 0.05 && a1 <= 0.95) ? 1.0 : 0.0;
            return {H * c.dV, std::sqrt(AMREX_D_TERM(gx*gx, + gy*gy, + gz*gz)) * c.dV};
        });

        const Real dx_f = FinestGeom(o).CellSize(0);
        columns.push_back({"N_I", (t[1] > 0.0) ? t[0] / t[1] / dx_f : 0.0});
    }

    WriteRow(o, "advection.csv", columns);
}

#else

// #######################################
// N-phase models (5Eq-N, 6Eq-N): four shapes of four materials
// #######################################

static_assert(NPHASE == 5, "This case sets five materials; build it with NPHASE = 5.");

// =======================================
// Initial condition
//
// Four shapes in material 1, each a tanh profile of width eta about a signed distance psi
// (positive inside): a slotted disk of material 2 (Zalesak, JCP 31, 1979), a square of
// material 3, a circle of material 4 and an equilateral triangle, apex down, of material 5.
// Material 1 is the complement, psi_1 = -max(psi_2, ..., psi_5), and the volume fractions
// are normalized to sum to one. Material k has density k. The pressure and velocity are
// uniform, and every material starts at the same pressure, so the interfaces only translate
// (Abgrall, JCP 125, 1996).
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

    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps = 1.0e-12;
    auto profile = [&](amrex::Real psi)
    {
        return 0.5*(1.0 + std::tanh( psi/(std::sqrt(2.0)*eta) )*(1.0 - eps) );
    };

    // Signed distance and density of each material
    amrex::Array<amrex::Real, NPHASE> Psi_phase, Alpha_phase, rho_phase;
    amrex::Real xc, yc, zc, r, L, X, Y;

    // ---------------------------------------
    // Material 2, slotted disk: radius 0.4, slot of width 0.125 cut up from the bottom to 0.65 deep
    // ---------------------------------------
    xc = 0.4;
    yc = -0.4;
    r  = 0.4;
    const amrex::Real slot_w = 0.125;
    const amrex::Real slot_d = 0.65;
    X = x - xc;
    Y = y - yc;
    const amrex::Real psi_disk   = r - std::sqrt(X*X + Y*Y);
    const amrex::Real psi_slot_x = (0.5 * slot_w) - std::abs(X);
    const amrex::Real psi_slot_y = (-r + slot_d) - Y;
    const amrex::Real psi_slot   = std::min(psi_slot_x, psi_slot_y);   // positive inside the slot
    Psi_phase[1] = std::min(psi_disk, -psi_slot);
    rho_phase[1] = 2.0;

    // ---------------------------------------
    // Material 3, square of side 0.4
    // ---------------------------------------
    xc = 0.5;
    yc = 0.5;
    L  = 0.4;
    Psi_phase[2] = 0.5 * L - std::max(std::abs(x - xc), std::abs(y - yc));
    rho_phase[2] = 3.0;

    // ---------------------------------------
    // Material 4, circle of radius 0.2 (a sphere in 3D)
    // ---------------------------------------
    xc = -0.5;
    yc = -0.5;
    zc =  0.0;
    r  =  0.2;
    Psi_phase[3] = r - std::sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
    rho_phase[3] = 4.0;

    // ---------------------------------------
    // Material 5, equilateral triangle, apex down: d1, d2 and d3 are the signed distances to
    // its three edges, which gives a height 2h/3 and a side 2L/3
    // ---------------------------------------
    xc = -0.5;
    yc =  0.5;
    L  = 0.75;
    const amrex::Real h = std::sqrt(3.0) * 0.5 * L;
    X = x - xc;
    Y = y - yc;
    const amrex::Real d1 =  (h/3.0) - Y;
    const amrex::Real d2 =  (std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    const amrex::Real d3 = -(std::sqrt(3.0)/2.0)*X + 0.5*Y + h/6.0;
    Psi_phase[4] = std::min(d1, std::min(d2, d3));
    rho_phase[4] = 5.0;

    // ---------------------------------------
    // Material 1 fills the rest
    // ---------------------------------------
    Psi_phase[0] = -std::max({Psi_phase[1], Psi_phase[2], Psi_phase[3], Psi_phase[4]});
    rho_phase[0] = 1.0;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_phase[i_Phase] = profile(Psi_phase[i_Phase]);
    }
    amrex::Real Sum_Alpha = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Sum_Alpha += Alpha_phase[i_Phase];
    }
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_phase[i_Phase] /= Sum_Alpha;
    }

    // ---------------------------------------
    // Uniform pressure and velocity
    // ---------------------------------------
    const amrex::Real u = Uniform::u;
    const amrex::Real v = Uniform::v;
    const amrex::Real w = 0.0;
    const amrex::Real P = Uniform::P;

#if (PHYSICS == FIVEEQS_NPHASE)
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha_phase[i_Phase]*rho_phase[i_Phase];
    }
#else
    amrex::Array<amrex::Real, NPHASE> Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha_phase[i_Phase]*rho_phase[i_Phase];
        Alpha_P  [i_Phase] = Alpha_phase[i_Phase]*P;
    }
#endif

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
#if (PHYSICS == FIVEEQS_NPHASE)
    State2Primitive(Primitive,
                    Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha_phase);
#else
    State2Primitive(Primitive,
                    Alpha_phase, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    P);
#endif
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/advection.csv with
//   V_1 ... V_5           volume (area in 2D) of each material, int alpha_k dV
//   alpha_min, alpha_max  extrema of the volume fractions over all materials
//   sum_err               max |1 - sum_k alpha_k|
//   dP                    max |P - P0|, with the mixture pressure P = sum_k alpha_k P_k in
//                         SIXEQS_IE_NPHASE
//   du                    max |u - u0| and |v - v0|
// where P0, u0 and v0 are the uniform state of the initial condition, which exact advection
// keeps. The integrals and the extrema are taken over the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const int ia = INDEX_VolumeFraction1;

    const auto V = Measure<NPHASE>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<NPHASE>
    {
        Values<NPHASE> v;
        for (int n = 0; n < NPHASE; ++n) { v[n] = U(c.i,c.j,c.k,ia+n) * c.dV; }
        return v;
    });

    const auto r = Measure(o, {Min, Max, Max, Max, Max}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        Real lo = 1.0e30, hi = -1.0e30, s = 0.0;
        for (int n = 0; n < NPHASE; ++n) {
            const Real a = U(c.i,c.j,c.k,ia+n);
            lo = std::min(lo, a);
            hi = std::max(hi, a);
            s += a;
        }

        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Conservative2Primitive(Prim, Cons, parm);
#if (PHYSICS == FIVEEQS_NPHASE)
        amrex::Array<Real, NPHASE> Alpha_rho, Alpha;
        Real P;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(Alpha_rho, AMREX_D_DECL(u, v, w), P, Alpha, Prim);
#else
        amrex::Array<Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
        Real P;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(Alpha, Alpha_rho, Alpha_P, AMREX_D_DECL(u, v, w), P, Prim);
#endif

        return {lo, hi, std::abs(1.0 - s), std::abs(P - Uniform::P),
                std::max(std::abs(u - Uniform::u), std::abs(v - Uniform::v))};
    });

    std::vector<std::pair<std::string, Real>> columns;
    for (int n = 0; n < NPHASE; ++n) { columns.push_back({"V_" + std::to_string(n+1), V[n]}); }
    columns.insert(columns.end(), {{"alpha_min", r[0]}, {"alpha_max", r[1]}, {"sum_err", r[2]},
                                   {"dP", r[3]}, {"du", r[4]}});
    WriteRow(o, "advection.csv", columns);
}

#endif

#endif  // IC_H_
