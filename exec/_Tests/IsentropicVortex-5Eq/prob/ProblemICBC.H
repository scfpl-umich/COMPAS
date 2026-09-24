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
// Isentropic vortex (Shu, Lecture Notes in Mathematics 1697, 1998)
//
// A vortex of strength epsilon = 5 centered at the origin in a uniform flow with
// rho = p = T = 1 and velocity (1, 1):
//   (du, dv) = epsilon/(2 pi) exp((1 - r^2)/2) (-y, x)
//   dT       = -(gamma - 1) epsilon^2/(8 gamma pi^2) exp(1 - r^2)
//   rho      = (1 + dT)^(1/(gamma - 1)),  p = rho^gamma
// Single fluid (alpha_1 = 1), ideal gas with gamma = prob.EOS.gamma_1.
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void isentropic_vortex(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       amrex::Real x,
                       amrex::Real y,
                       Parm const& parm)
{
    const amrex::Real epsilon = 5.0;           // vortex strength
    const amrex::Real gamma   = parm.gamma_1;
    const amrex::Real u_inf   = 1.0;           // uniform flow at 45 degrees
    const amrex::Real v_inf   = 1.0;

    const amrex::Real r2   = x*x + y*y;
    const amrex::Real exp1 = std::exp(0.5*(1.0 - r2));
    const amrex::Real exp2 = std::exp(1.0 - r2);

    const amrex::Real dvx = - epsilon/(2.0*PI) * exp1 * y;
    const amrex::Real dvy =   epsilon/(2.0*PI) * exp1 * x;
    const amrex::Real dT  = - (gamma - 1.0) * epsilon*epsilon / (8.0 * gamma * PI * PI) * exp2;

    const amrex::Real rho = std::pow(1.0 + dT, 1.0/(gamma - 1.0));
    const amrex::Real u   = u_inf + dvx;
    const amrex::Real v   = v_inf + dvy;
    const amrex::Real w   = 0.0;
    const amrex::Real p   = std::pow(rho, gamma);

    // Convert to primitive and then conservative
    const amrex::Real Alpha1 = 1.0;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1 * rho,
                    (1.0 - Alpha1) * rho,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha1);
    Primitive2Conservative(Conservative, Primitive, parm);
}


// =======================================
// Initial condition: the vortex at the origin
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
    amrex::ignore_unused(z, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    isentropic_vortex(Conservative, x, y, parm);
}


// =======================================
// Exact solution, used with CONVERGENCE for the error norms at the end of the run
//
// The initial vortex translated by (t, t) and wrapped into the periodic box
// [-5, 5]^2 (geometry.prob_lo and prob_hi in prob/inputs), so that it returns to the
// initial condition after one period, t = 10.
// =======================================

#undef  CONVERGENCE          // Macros.H defines it as false by default
#define CONVERGENCE true
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void exact_solution(amrex::Array<amrex::Real, NSTATE>& Conservative,
                    amrex::Real x,
                    amrex::Real y,
                    amrex::Real z,
                    amrex::Real t,
                    AMREX_D_DECL(amrex::Real dx,
                                 amrex::Real dy,
                                 amrex::Real dz),
                    AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                 amrex::Real dy_FinestLevel,
                                 amrex::Real dz_FinestLevel),
                    Parm const& parm)
{
    amrex::ignore_unused(z, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    // Periodic image of the point, traced back to t = 0 along the uniform flow (1, 1)
    const amrex::Real L  = 10.0;
    const amrex::Real x0 = x - t;
    const amrex::Real y0 = y - t;
    const amrex::Real xs = x0 - L * std::floor((x0 + 0.5*L) / L);
    const amrex::Real ys = y0 - L * std::floor((y0 + 0.5*L) / L);

    isentropic_vortex(Conservative, xs, ys, parm);
}

#endif  // IC_H_
