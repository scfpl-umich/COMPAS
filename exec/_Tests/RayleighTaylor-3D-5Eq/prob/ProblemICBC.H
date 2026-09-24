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
// Three-dimensional single-mode Rayleigh-Taylor instability. Material 1 (heavy) fills
// y > a0 cos(2 pi x/lambda) cos(2 pi z/lambda), with lambda = 1, and material 2 (light) is
// below it. Both are at rest. The interface is a tanh profile of width
// eta = PhaseField.Eta_Multiplier finest cells, and the pressure is hydrostatic in each
// column, P = P0 + g_y int_0^y rho dy', so it is continuous across the interface.
// =======================================

// log cosh(s), without overflow for large |s|
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real LogCosh(amrex::Real s)
{
    s = std::abs(s);
    return s + std::log1p(std::exp(-2.0 * s)) - std::log(2.0);
}

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
    // Interface: tanh profile about y = a0 cos(2 pi x/lambda) cos(2 pi z/lambda)
    // ---------------------------------------
    const amrex::Real eta    = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps    = 1.0e-12;
    const amrex::Real lambda = 1.0;
    const amrex::Real h      = parm.a0 * std::cos(2.0 * PI * x / lambda) * std::cos(2.0 * PI * z / lambda);
    const amrex::Real Psi    = y - h;

    const amrex::Real delta  = std::sqrt(2.0) * eta;
    const amrex::Real Alpha1 = 0.5 * (1.0 + std::tanh(Psi / delta) * (1.0 - eps));
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Fluids at rest in hydrostatic balance along the column, dP/dy = rho g_y with
    // rho = rho2 + (rho1 - rho2) alpha_1, integrated exactly from P = P0 at y = 0 using
    // int_0^y alpha_1 dy' = y/2 + (1 - eps) (delta/2) [log cosh((y - h)/delta) - log cosh(-h/delta)]
    // ---------------------------------------
    const amrex::Real Alpha1_rho1 = Alpha1 * parm.rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * parm.rho2;

    const amrex::Real Int_Alpha1 = 0.5 * y + 0.5 * (1.0 - eps) * delta * (LogCosh(Psi / delta) - LogCosh(-h / delta));

    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    const amrex::Real P = parm.P0 + parm.gy * (parm.rho2 * y + (parm.rho1 - parm.rho2) * Int_Alpha1);

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
// Adds rho g to the momentum equations and rho g . u to the energy equation. The mass and
// volume-fraction equations have no source.
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
// With run.user_output_int = N, COMPAS calls this every N coarse steps, and otherwise with
// each plotfile. Each call appends one line to <plot directory>/out1D_<case name>.csv with
//   I_L   interface position, where alpha_1 crosses 0.5 in y, at x = z = 0 (the bubble tip)
//   I_R   the same at x = lambda/2, z = 0 (the spike tip)
//   A     amplitude, |I_L - I_R| / 2
// Crossings are taken on the finest level, in the plane z = 0.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // One ghost cell, for the -y neighbor in the crossing
    const auto v = Measure<2>(o, Max, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<2>
    {
        // alpha_1 = 0.5 crossing toward the -y neighbor
        const Real aL  = U(c.i, c.j-1, c.k, INDEX_VolumeFraction1);
        const Real aR  = U(c.i, c.j,   c.k, INDEX_VolumeFraction1);
        const int  klo = (AMREX_SPACEDIM == 3) ? c.domain.smallEnd(AMREX_SPACEDIM-1) : 0;
        if (c.k != klo || !(aL <= 0.5 && aR >= 0.5 && aR > aL)) { return {-big, -big}; }
        const Real cross = c.x[1] - c.dx[1] + c.dx[1] * (0.5 - aL) / (aR - aL);
        return {(c.i == c.domain.smallEnd(0)) ? cross : -big, (c.i == c.domain.bigEnd(0)) ? cross : -big};
    }, FinestLevel);

    WriteRow(o, "out1D_" + o.case_name + ".csv", {{"I_L", v[0]}, {"I_R", v[1]}, {"A", 0.5 * std::abs(v[0] - v[1])}});
}

#endif  // IC_H_
