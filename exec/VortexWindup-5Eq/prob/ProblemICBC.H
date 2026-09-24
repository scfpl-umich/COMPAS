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
// A steady isentropic vortex at rest (Shu, Lecture Notes in Mathematics 1697, 1998) winds up a
// planar interface between two materials. With r' = |x - c|/rc and e = exp((1 - r'^2)/2),
//   (u, v) = beta/(2 pi) e (-(y - yc), x - xc)/rc
//   T/T_inf = 1 - (gamma - 1) beta^2/(8 gamma pi^2 T_inf) e^2,   T_inf = p0/rho0
//   rho = rho0 (T/T_inf)^(1/(gamma - 1)),  p = p0 (T/T_inf)^(gamma/(gamma - 1))
// The pressure gradient balances the centripetal acceleration, so the flow is steady. Both
// materials are the same ideal gas at the same density, so the interface is passive: each point
// turns about the center at the angular velocity
//   Omega(r) = beta/(2 pi rc) exp((1 - r'^2)/2)
// and the exact volume fraction at time t is the initial one turned by -Omega(r) t at each radius.
// Omega falls off with r, so the interface winds up into a double spiral whose arms thin in time.
//
// Material 2 fills the half plane n.(x - c) > offset, or with prob.width > 0 the slab
// 0 < n.(x - c) - offset < width, n = (cos angle, sin angle). Each cell holds the exact cell
// average of alpha_2 (the area of the cell on the side of material 2, by polygon clipping), bounded
// to [alpha_min, 1 - alpha_min]; the flow is the average of the vortex over 3 x 3 Gauss points.
// Set run.QuadratureIC = 0, since the cell averages are computed here.
// =======================================

namespace VortexWindup
{
using amrex::Real;

// Convex polygon, clipped by half planes a x + b y <= c
struct Poly
{
    int  n = 0;
    Real x[8], y[8];
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Clip (Poly& P, Real a, Real b, Real c)
{
    Poly Q;
    for (int m = 0; m < P.n; ++m) {
        const int  q  = (m + 1) % P.n;
        const Real fm = a * P.x[m] + b * P.y[m] - c;
        const Real fq = a * P.x[q] + b * P.y[q] - c;
        if (fm <= 0.0) { Q.x[Q.n] = P.x[m]; Q.y[Q.n] = P.y[m]; Q.n++; }
        if ((fm < 0.0 && fq > 0.0) || (fm > 0.0 && fq < 0.0)) {
            const Real s = fm / (fm - fq);
            Q.x[Q.n] = P.x[m] + s * (P.x[q] - P.x[m]);
            Q.y[Q.n] = P.y[m] + s * (P.y[q] - P.y[m]);
            Q.n++;
        }
    }
    P = Q;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real Area (Poly const& P)
{
    Real A = 0.0;
    for (int m = 0; m < P.n; ++m) {
        const int q = (m + 1) % P.n;
        A += P.x[m] * P.y[q] - P.x[q] * P.y[m];
    }
    return 0.5 * std::abs(A);
}

// Exact cell average of the initial alpha_2 over the cell of center (x, y) and size dx x dy
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real CellAlpha2 (Real x, Real y, Real dx, Real dy, Parm const& parm)
{
    Poly P;
    P.n = 4;
    P.x[0] = x - 0.5*dx; P.y[0] = y - 0.5*dy;
    P.x[1] = x + 0.5*dx; P.y[1] = y - 0.5*dy;
    P.x[2] = x + 0.5*dx; P.y[2] = y + 0.5*dy;
    P.x[3] = x - 0.5*dx; P.y[3] = y + 0.5*dy;
    const Real nc = parm.nx * parm.xc + parm.ny * parm.yc;
    // n.x >= nc + offset
    Clip(P, -parm.nx, -parm.ny, -(nc + parm.offset));
    // n.x <= nc + offset + width
    if (parm.width > 0.0) { Clip(P, parm.nx, parm.ny, nc + parm.offset + parm.width); }
    return (P.n > 2) ? Area(P) / (dx * dy) : 0.0;
}

// Initial alpha_2 at a point: 1 on the side of material 2, 0 elsewhere
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real PointAlpha2 (Real x, Real y, Parm const& parm)
{
    const Real d = parm.nx * (x - parm.xc) + parm.ny * (y - parm.yc) - parm.offset;
    return (d > 0.0 && (parm.width <= 0.0 || d < parm.width)) ? 1.0 : 0.0;
}

// Exact alpha_2 at time t, averaged over n x n points of the cell: each point is turned back by
// the angle Omega(r) t about the center
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real ExactAlpha2 (Real x, Real y, Real dx, Real dy, Real t, int n, Parm const& parm)
{
    Real a = 0.0;
    for (int jj = 0; jj < n; ++jj) {
        for (int ii = 0; ii < n; ++ii) {
            const Real X  = x + ((ii + 0.5) / n - 0.5) * dx - parm.xc;
            const Real Y  = y + ((jj + 0.5) / n - 0.5) * dy - parm.yc;
            const Real r2 = (X*X + Y*Y) / (parm.rc * parm.rc);
            const Real th = -parm.beta / (2.0 * M_PI * parm.rc) * std::exp(0.5 * (1.0 - r2)) * t;
            const Real c  = std::cos(th), s = std::sin(th);
            a += PointAlpha2(parm.xc + c*X - s*Y, parm.yc + s*X + c*Y, parm);
        }
    }
    return a / Real(n * n);
}

// Density, velocity and pressure of the vortex at the point (x, y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Vortex (Real x, Real y, Parm const& parm, Real& rho, Real& u, Real& v, Real& p)
{
    const Real g    = parm.gamma_1;
    const Real X    = x - parm.xc, Y = y - parm.yc;
    const Real r2   = (X*X + Y*Y) / (parm.rc * parm.rc);
    const Real e    = std::exp(0.5 * (1.0 - r2));
    const Real Tinf = parm.p0 / parm.rho0;
    const Real T    = 1.0 - (g - 1.0) * parm.beta * parm.beta / (8.0 * g * M_PI * M_PI * Tinf) * e * e;
    rho = parm.rho0 * std::pow(T, 1.0 / (g - 1.0));
    p   = parm.p0   * std::pow(T, g / (g - 1.0));
    u   = -parm.beta / (2.0 * M_PI) * e * Y / parm.rc;
    v   =  parm.beta / (2.0 * M_PI) * e * X / parm.rc;
}

}  // namespace VortexWindup


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
    amrex::ignore_unused(z, AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));
#if (AMREX_SPACEDIM == 3)
    amrex::ignore_unused(dz);
#endif

    // ---------------------------------------
    // Volume fractions: the exact cell average of the planar interface
    // ---------------------------------------
    const amrex::Real a2     = VortexWindup::CellAlpha2(x, y, dx, dy, parm);
    const amrex::Real Alpha2 = parm.alpha_min + (1.0 - 2.0 * parm.alpha_min) * a2;
    const amrex::Real Alpha1 = 1.0 - Alpha2;

    // ---------------------------------------
    // Vortex: average of the conservative state over 3 x 3 Gauss points. The materials are the
    // same gas at the same density, so alpha_k rho_k = Alpha_k rho at each point.
    // ---------------------------------------
    const amrex::Real q[3] = {-0.5 * std::sqrt(0.6), 0.0, 0.5 * std::sqrt(0.6)};
    const amrex::Real wq[3] = {5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};
    for (int n = 0; n < NSTATE; ++n) { Conservative[n] = 0.0; }
    for (int jq = 0; jq < 3; ++jq) {
        for (int iq = 0; iq < 3; ++iq) {
            amrex::Real rho, u, v, P;
            VortexWindup::Vortex(x + q[iq] * dx, y + q[jq] * dy, parm, rho, u, v, P);
            const amrex::Real w = 0.0;
            amrex::ignore_unused(w);

            // ---------------------------------------
            // Convert to primitive and then conservative
            // ---------------------------------------
            amrex::Array<amrex::Real, NSTATE> Primitive, Cq;
            State2Primitive(Primitive,
                            Alpha1 * rho, Alpha2 * rho,
               AMREX_D_DECL(u, v, w),
                            P,
                            Alpha1);
            Primitive2Conservative(Cq,
                                   Primitive,
                                   parm);
            for (int n = 0; n < NSTATE; ++n) { Conservative[n] += wq[iq] * wq[jq] * Cq[n]; }
        }
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends one
// line to <plot directory>/windup.csv with
//   V2              volume of material 2, int alpha_2 dA (conserved)
//   a1_min, a1_max  extremes of alpha_1
//   A_mix           area of the mixed cells, 0.01 < alpha_1 < 0.99
//   L1              int |alpha_2 - alpha_2,exact| dA, with the exact alpha_2 averaged over
//                   prob.nsub x prob.nsub points per cell (analysis/windup.py computes it on
//                   finer samples from the plotfiles)
// The integrals are taken over the composite grid, i.e. each cell not covered by a finer level.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real t = Time(o);

    const auto m = Measure(o, {Sum, Min, Max, Sum, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        const Real a1 = U(c.i, c.j, c.k, INDEX_VolumeFraction1);
        const Real a2 = 1.0 - a1;
        const Real ae = VortexWindup::ExactAlpha2(c.x[0], c.x[1], c.dx[0], c.dx[1], t, parm.nsub, parm);
        return {a2 * c.dV, a1, a1, (a1 > 0.01 && a1 < 0.99) ? c.dV : 0.0, std::abs(a2 - ae) * c.dV};
    });

    WriteRow(o, "windup.csv", {{"V2", m[0]}, {"a1_min", m[1]}, {"a1_max", m[2]}, {"A_mix", m[3]}, {"L1", m[4]}});
}

#endif  // IC_H_
