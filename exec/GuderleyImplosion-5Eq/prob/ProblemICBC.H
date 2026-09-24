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
// Guderley implosion onto a perturbed interface, on the quarter plane x, y >= 0.
// A strong cylindrical shock of Mach number prob.Mach is at r = prob.Rs, with the Guderley
// similarity solution behind it (see Guderley.H), and converges on a heavier inner gas
// (phase 2) whose interface is R(theta) = prob.R0 + prob.a0 cos(prob.mode theta). Both
// gases ahead of the shock are at rest at the pressure prob.P0.
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
    // Gases at rest ahead of the shock
    // ---------------------------------------
    amrex::Real rho1 = parm.rho1;
    amrex::Real rho2 = parm.rho2;
    amrex::Real u_r  = 0.0;
    amrex::Real P    = parm.P0;

    // ---------------------------------------
    // Behind the shock, the similarity solution at t = 0:
    //   u = -(r/tau_c) U(xi), c^2 = (r/tau_c)^2 C(xi), rho = rho1 G(xi), xi = r/Rs
    // ---------------------------------------
    if (r > parm.Rs){
        Guderley::State(parm.ss, r, 0.0, rho1, u_r, P);
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
// Optional: boundary conditions
//
// Behind a strong converging shock the gas flows inward supersonically, so every wave at
// the outer boundary comes from outside, and the boundary needs the exact inflow. With
// run.hi_bc = 1 (inflow), the ghost cells there are filled with the similarity solution at
// the current time. This is valid until the shock reflected from the center reaches the
// boundary.
// =======================================

#define USER_BOUNDARY_FUNC
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_User(const IntVect& iv,
                       Array4<Real> const& q,
                       const int dcomp,
                       const int numcomp,
                       GeometryData const& geom,
                       const Real time,
                       const BCRec* bcr,
                       const int bcomp,
                       const int orig_comp,
                       Parm const* parm)
{
    const auto hi  = amrex::ubound(geom.domain);
    const bool xhi = (iv[0] > hi.x) && (bcr[bcomp].hi(0) == BCType::ext_dir);
    const bool yhi = (iv[1] > hi.y) && (bcr[bcomp].hi(1) == BCType::ext_dir);
    if (!xhi && !yhi) { return; }

    // Fields other than the whole state (Phase-Field and level-set buffers) get a
    // zero-gradient copy of the nearest cell inside the domain
    if (dcomp != 0 || orig_comp != 0 || numcomp < NSTATE) {
        IntVect ib = iv;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ib[d] = amrex::max(geom.domain.smallEnd(d), amrex::min(iv[d], geom.domain.bigEnd(d)));
        }
        for (int n = dcomp; n < dcomp + numcomp; ++n) { q(iv, n) = q(ib, n); }
        return;
    }

    const Real x = geom.ProbLo(0) + (iv[0] + 0.5) * geom.CellSize(0);
    const Real y = geom.ProbLo(1) + (iv[1] + 0.5) * geom.CellSize(1);
    const Real r = std::sqrt(x*x + y*y);

    Real rho = parm->rho1, u_r = 0.0, P = parm->P0;
    Guderley::State(parm->ss, r, time, rho, u_r, P);

    // Pure air, with the same trace of phase 2 as the initial condition far from the interface
    const Real Alpha2 = 0.5e-8;
    Array<Real, NSTATE> Primitive, Conservative;
    State2Primitive(Primitive, (1.0 - Alpha2) * rho, Alpha2 * parm->rho2,
       AMREX_D_DECL(u_r * x / r, u_r * y / r, 0.0), P, 1.0 - Alpha2);
    Primitive2Conservative(Conservative, Primitive, *parm);
    for (int n = dcomp; n < dcomp + numcomp; ++n) {
        const int s = orig_comp + n - dcomp;
        if (s < NSTATE) { q(iv, n) = Conservative[s]; }
    }
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
//   r_shock  smallest radius where P > 10 prob.P0, i.e. the converging shock until it focuses
// The integrals are taken over the composite grid, i.e. each cell not covered by a finer level.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    Parm const& hparm   = HostParm(o);
    const int   mode    = hparm.mode;
    const Real  P_shock = 10.0 * hparm.P0;
    const Real  Theta   = 0.5 * PI;   // quarter plane
    const Real  big     = 1.0e30;

    // One ghost cell, for the alpha_1 = 0.5 crossings between neighbors
    const auto v = Measure(o, {Sum, Sum, Min, Max, Min}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real r      = std::sqrt(c.x[0]*c.x[0] + c.x[1]*c.x[1]);
        const Real theta  = std::atan2(c.x[1], c.x[0]);
        const Real alpha1 = U(i,j,k,INDEX_VolumeFraction1);
        const Real alpha2 = 1.0 - alpha1;

        // Interface: alpha_1 - 0.5 changes sign toward the +x or +y neighbor on the finest level
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
        AMREX_D_TERM(Real ux;, Real uy;, Real uz;)
        Primitive2State(a1r1, a2r2, AMREX_D_DECL(ux, uy, uz), P, a1, Prim);

        return {alpha2 * c.dV, alpha2 * std::cos(mode * theta) * c.dV, r_lo, r_hi, (P > P_shock) ? r : big};
    });

    // The similarity solution, for comparing with theory in Python
    if (Step(o) == 0) { Guderley::Write(hparm.ss, hparm.P0, Directory(o)); }

    const Real R_eq = std::sqrt(4.0 * v[0] / PI);
    WriteRow(o, "interface.csv", {{"R_eq", R_eq}, {"a_" + std::to_string(mode), 2.0 * v[1] / (R_eq * Theta)},
                                  {"h", 0.5 * (v[3] - v[2])}, {"r_shock", (v[4] < big) ? v[4] : 0.0}});
}

#endif  // IC_H_
