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


// Uniform pressure and velocity of the flow, shared by the initial condition and the diagnostics
namespace Uniform {
    constexpr amrex::Real P = 1.0;
    constexpr amrex::Real u = 1.0;
    constexpr amrex::Real v = 1.0;
}


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

    const amrex::Real u = Uniform::u;
    const amrex::Real v = Uniform::v;
    const amrex::Real w = 0.0;
    const amrex::Real P = Uniform::P;

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
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/advection.csv with
//   V_2                   volume (area in 2D) of material 2, int (1 - alpha_1) dV
//   alpha_min, alpha_max  extrema of alpha_1
//   dP                    max |P - P0|
//   du                    max |u - u0| and |v - v0|
// where P0, u0 and v0 are the uniform state of the initial condition, which exact advection
// keeps. The integral and the extrema are taken over the composite grid.
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
        Real a1r1, a2r2, P, a1;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(a1r1, a2r2, AMREX_D_DECL(u, v, w), P, a1, Prim);

        return {(1.0 - a1) * c.dV, a1, a1, std::abs(P - Uniform::P),
                std::max(std::abs(u - Uniform::u), std::abs(v - Uniform::v))};
    });

    WriteRow(o, "advection.csv", {{"V_2", r[0]}, {"alpha_min", r[1]}, {"alpha_max", r[2]},
                                  {"dP", r[3]}, {"du", r[4]}});
}

#endif  // IC_H_
