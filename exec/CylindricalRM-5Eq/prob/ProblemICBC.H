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
// Cylindrical Richtmyer-Meshkov problem on the quarter plane x, y >= 0.
// A cylindrical shock of Mach number prob.Mach starts at r = prob.Rs in the outer
// gas (phase 1) and converges on a heavier inner gas (phase 2) whose interface is
//     R(theta) = prob.R0 + prob.a0 cos(prob.mode theta).
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
    // Polar coordinates and interface width
    // ---------------------------------------
    const amrex::Real r     = std::sqrt(x*x + y*y);
    const amrex::Real theta = std::atan2(y, x);
    const amrex::Real eta   = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real eps   = 1.0e-8;

    // ---------------------------------------
    // Perturbed interface: phase 2 inside R(theta)
    // ---------------------------------------
    const amrex::Real R_Interface = parm.R0 + parm.a0 * std::cos(parm.mode * theta);
    const amrex::Real Psi         = R_Interface - r;   // > 0 inside the interface
    const amrex::Real Alpha2      = 0.5 * (1.0 + std::tanh(Psi / (std::sqrt(2.0) * eta)) * (1.0 - eps));
    const amrex::Real Alpha1      = 1.0 - Alpha2;

    // ---------------------------------------
    // Unshocked gases at rest, P = 1
    // ---------------------------------------
    amrex::Real rho1 = parm.rho1;
    amrex::Real rho2 = parm.rho2;
    amrex::Real u_r  = 0.0;
    amrex::Real P    = 1.0;

    // ---------------------------------------
    // Converging shock: post-shock state of phase 1 outside r = Rs,
    // from the Rankine-Hugoniot relations, moving inward
    // ---------------------------------------
    if (r > parm.Rs){
        amrex::Real shock_param[5];
        StiffenedShockConditions(parm.Mach, 0.0, P, parm.rho1, parm.gamma_1, parm.pinf_1, shock_param);
        rho1 = shock_param[0];
        u_r  = -shock_param[1];
        P    = shock_param[2];
    }
    const amrex::Real u = (r > 0.0) ? u_r * x / r : 0.0;
    const amrex::Real v = (r > 0.0) ? u_r * y / r : 0.0;
    const amrex::Real w = 0.0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1 * rho1, Alpha2 * rho2,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Optional: in-situ diagnostics
//
// With run.user_output_int = N > 0 in the inputs file, COMPAS calls this function every
// N coarse steps, which records the interface far more often than the plotfiles do.
// Each call appends one line to <plot directory>/interface.csv with
//   R_eq     equivalent radius of the inner gas, from its area A = pi R_eq^2 / 4
//   a_n      signed amplitude of mode n, 2 int alpha_2 cos(n theta) dA / (R_eq Theta)
//   h        half the radial extent of the alpha_1 = 0.5 contour on the finest level
//   r_shock  smallest radius where P > 1.5, i.e. the converging shock until it focuses
// The integrals are over the composite grid, the cells that no finer level covers, and
// use the helpers of include/Tools/UserOutput.H.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const int  mode    = HostParm(o).mode;
    const Real Theta   = 0.5 * PI;   // quarter plane
    const Real big     = 1.0e30;
    const Real P_shock = 1.5;        // between the pre- and post-shock pressures

    // One pass over the composite grid, with one ghost cell so that each cell sees its
    // neighbors: two sums, the extent of the interface in r, and the shock radius
    const auto m = Measure(o, {Sum, Sum, Min, Max, Min}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real r      = std::sqrt(c.x[0]*c.x[0] + c.x[1]*c.x[1]);
        const Real theta  = std::atan2(c.x[1], c.x[0]);
        const Real alpha1 = U(i,j,k,INDEX_VolumeFraction1);
        const Real alpha2 = 1.0 - alpha1;

        // Interface on the finest level: alpha_1 - 0.5 changes sign toward +x or +y
        Real r_lo = big, r_hi = 0.0;
        if (c.finest) {
            const Real s  = alpha1 - 0.5;
            const Real sx = U(i+1,j,k,INDEX_VolumeFraction1) - 0.5;
            const Real sy = U(i,j+1,k,INDEX_VolumeFraction1) - 0.5;
            if (s * sx <= 0.0 || s * sy <= 0.0) { r_lo = r; r_hi = r; }
        }

        // Pressure, for the shock position
        amrex::Array<Real, NSTATE> Cons, Prim;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(i,j,k,n); }
        Conservative2Primitive(Prim, Cons, parm);
        Real a1r1, a2r2, a1, P;
        AMREX_D_TERM(Real u;, Real v;, Real w;)
        Primitive2State(a1r1, a2r2, AMREX_D_DECL(u, v, w), P, a1, Prim);

        return {alpha2 * c.dV, alpha2 * std::cos(mode * theta) * c.dV,
                r_lo, r_hi, (P > P_shock) ? r : big};
    });

    const Real R_eq    = std::sqrt(4.0 * m[0] / PI);
    const Real a_n     = 2.0 * m[1] / (R_eq * Theta);
    const Real h       = 0.5 * (m[3] - m[2]);
    const Real r_shock = (m[4] < big) ? m[4] : 0.0;

    WriteRow(o, "interface.csv", {{"R_eq", R_eq}, {"a_" + std::to_string(mode), a_n},
                                  {"h", h}, {"r_shock", r_shock}});
}

#endif  // IC_H_
