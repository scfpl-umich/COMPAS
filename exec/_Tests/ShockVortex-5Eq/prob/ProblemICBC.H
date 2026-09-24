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
// Initial condition: shock-vortex interaction (Rault, Chiavassa and Donat,
// J. Sci. Comput. 19, 2003)
//
// A stationary Mach 1.5 normal shock at x = 0.5, with rho = p = 1 (and R = 1) upstream.
// A composite vortex of Mach number 0.9 centered at (0.25, 0.5) is superimposed on the
// upstream flow. Its azimuthal velocity grows linearly up to r = a, then decays as
// (r - b^2/r) to zero at r = b, and the temperature follows from radial equilibrium with
// isentropic pressure. Single fluid (material 2), ideal gas with gamma = prob.EOS.gamma_2.
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

    const Real gamma = parm.gamma_2;      // the fluid is material 2
    const Real gm1   = gamma - 1.0;
    const Real gp1   = gamma + 1.0;
    const Real R     = 1.0;

    // ---------------------------------------
    // Upstream state, and the downstream state from the normal-shock relations
    // ---------------------------------------
    const Real M_s   = 1.5;
    const Real Ms2   = M_s * M_s;
    const Real rho_u = 1.0;
    const Real p_u   = 1.0;
    const Real u_u   = M_s * std::sqrt(gamma * p_u / rho_u);
    const Real v_u   = 1.0e-20;
    const Real w_u   = 1.0e-20;
    const Real t_u   = p_u / (rho_u * R);

    const Real rho_d = rho_u * gp1 * Ms2 / (2.0 + gm1 * Ms2);
    const Real u_d   = rho_u * u_u / rho_d;
    const Real p_d   = p_u * (1.0 + (2.0 * gamma / gp1) * (Ms2 - 1.0));
    const Real t_d   = p_d / (rho_d * R);

    Real P, T, u, v, w;
    if (x <= 0.5) {
        P = p_u;  T = t_u;
        u = u_u;  v = v_u;  w = w_u;
    } else {
        P = p_d;  T = t_d;
        u = u_d;  v = v_u;  w = w_u;
    }

    // ---------------------------------------
    // Vortex: center (x_c, y_c), core radius a, outer radius b, peak speed v_m at r = a
    // ---------------------------------------
    const Real x_c = 0.25;
    const Real y_c = 0.5;
    const Real a   = 0.075;
    const Real b   = 0.175;

    const Real M_v   = 0.9;
    const Real v_m   = M_v * std::sqrt(gamma);
    const Real a2    = a * a;
    const Real b2    = b * b;
    const Real denom = (a2 - b2);
    const Real coeff = gm1 / (R * gamma);

    const Real dX = x - x_c;
    const Real dY = y - y_c;
    const Real r2 = dX*dX + dY*dY;
    const Real r  = std::sqrt(r2);

    if (r <= b) {
        const Real inv_r = 1.0 / (r + 1e-16);
        const Real sin_t = dY * inv_r;
        const Real cos_t = dX * inv_r;
        const Real vm2   = (v_m * a / denom) * (v_m * a / denom);

        if (r <= a) {
            // Solid-body core
            const Real mag = v_m * r / a;
            u -= mag * sin_t;
            v += mag * cos_t;

            // Temperature at r = a from the outer solution, then integrated inward
            const Real radial_a = -2.0*b2*std::log(b) - 0.5*a2 + 2.0*b2*std::log(a) + 0.5*b2*b2/a2;
            const Real t_a      = t_u - coeff * vm2 * radial_a;
            const Real radial_term = 0.5 * (1.0 - r2/a2);
            T = t_a - coeff * v_m * v_m * radial_term;
        } else {
            // Outer ring, with zero velocity at r = b
            const Real mag = v_m * a * (r - b2*inv_r) / denom;
            u -= mag * sin_t;
            v += mag * cos_t;

            const Real radial_term = -2.0*b2*std::log(b) - 0.5*r2 + 2.0*b2*std::log(r) + 0.5*b2*b2/r2;
            T = t_u - coeff * vm2 * radial_term;
        }

        // Isentropic pressure
        P = p_u * std::pow(T / t_u, gamma / gm1);
    }

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    const Real rho    = P / (R * T);
    const Real Alpha1 = 1.0e-12;

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1 * rho, (1.0 - Alpha1) * rho,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative, Primitive, parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/vortex.csv with
//   x_v, y_v   position of the vortex core, the cell of minimum pressure
//   p_min      the minimum pressure
// The minimum is taken over the uncovered cells of the composite grid.
// =======================================

// Pressure of cell (i,j,k)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real cell_pressure(amrex::Array4<amrex::Real const> const& Ua, int i, int j, int k, Parm const& parm)
{
    amrex::Array<amrex::Real, NSTATE> Cons, Prim;
    for (int n = 0; n < NSTATE; ++n) { Cons[n] = Ua(i,j,k,n); }
    Conservative2Primitive(Prim, Cons, parm);
    return Prim[2 + AMREX_SPACEDIM];    // pressure in the primitive vector
}

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // The minimum pressure, then the position of the cell where it is attained
    const Real p_min = Measure<1>(o, Min, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<1>
    {
        return {cell_pressure(U, c.i, c.j, c.k, parm)};
    })[0];
    const auto x_v = Measure<2>(o, Max, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<2>
    {
        if (cell_pressure(U, c.i, c.j, c.k, parm) > p_min) { return {-big, -big}; }
        return {c.x[0], c.x[1]};
    });
    WriteRow(o, "vortex.csv", {{"x_v", x_v[0]}, {"y_v", x_v[1]}, {"p_min", p_min}});
}

#endif  // IC_H_
