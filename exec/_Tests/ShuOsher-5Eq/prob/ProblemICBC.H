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
// Shu-Osher shock tube (Shu and Osher, J. Comput. Phys. 83, 1989).
// A Mach 3 shock at x = -4 runs into the density perturbation rho = 1 + 0.2 sin(5x).
// Behind the shock the perturbation is compressed into short entropy waves, which are
// sensitive to the dissipation of the reconstruction scheme. Single fluid, gamma = 1.4.
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
    amrex::ignore_unused(y, z, AMREX_D_DECL(dx, dy, dz));

    // ---------------------------------------
    // Post-shock state on the left, perturbed density at rest on the right
    // ---------------------------------------
    const bool left = (x < -4.0);
    const amrex::Real rho = left ? 3.857143 : 1.0 + 0.2 * std::sin(5.0 * x);
    const amrex::Real u   = left ? 2.629369 : 0.0;
    const amrex::Real v   = 0.0;
    const amrex::Real w   = 0.0;
    const amrex::Real P   = left ? 10.33333 : 1.0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    const amrex::Real Alpha1 = 1.0;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*rho, (1.0 - Alpha1)*rho,
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
// one line to <plot directory>/shock.csv with
//   x_shock    shock position, the furthest point in x where the pressure crosses 1.5 times
//              the pressure ahead of the shock (1), found on the finest level
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real P_thr = 1.5 * 1.0;

    // Finest level only, with one ghost cell for the crossing with the +x neighbor
    const auto v = Measure(o, {Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<1>
    {
        amrex::Array<Real, NSTATE> Uc, Ue;
        for (int n = 0; n < NSTATE; ++n) { Uc[n] = U(c.i,c.j,c.k,n); Ue[n] = U(c.i+1,c.j,c.k,n); }
        Real Pc, Pe;
        Conservative2Pressure(Pc, Uc, parm);
        Conservative2Pressure(Pe, Ue, parm);
        const bool cross = (Pc - P_thr) * (Pe - P_thr) <= 0.0 && Pc != Pe;
        return {cross ? c.x[0] + c.dx[0] * (P_thr - Pc) / (Pe - Pc) : -1.0e30};
    }, FinestLevel);
    WriteRow(o, "shock.csv", {{"x_shock", v[0]}});
}

#endif  // IC_H_
