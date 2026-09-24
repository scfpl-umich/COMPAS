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
// Pseudo-steady refraction of a planar shock in air (material 1) over a wedge of water
// (material 2), after Anbu Serene Raj et al. (J. Fluid Mech. 998, A49, 2024). The water fills
// y < (x - prob.x0) tan(theta_w) above the bottom wall, with its apex at (prob.x0, 0), and
// beta = 90 deg - theta_w (prob.beta) is the angle between the shock and the interface. Both
// fluids are at rest at the pressure prob.P0. A shock of Mach number prob.Mach travels in +x
// through the air from x = prob.x_s, with the Rankine-Hugoniot state behind it (Parm.H). The
// interface starts sharp: the water volume fraction of each cell is the fraction of its area
// that lies in the wedge (WaterFraction), bounded away from 0 and 1 by eps/2, so the initial
// data are exact cell averages. With run.QuadratureIC = 0 (prob/inputs) the function is called
// once per cell, at its center (x, y), with the cell size (dx, dy).
// =======================================

// Fraction of the cell of size dx by dy centered at (x, y) that lies in the water, where
// (x - x0) sin(theta_w) - y cos(theta_w) > 0, for 0 < theta_w < 90 deg. With a the value of
// this distance at the lower-right corner of the cell, the corner deepest in the water, the
// area is [F(a) - F(a - dx sin) - F(a - dy cos) + F(a - dx sin - dy cos)] / (2 sin cos),
// with F(s) = max(s, 0)^2.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real WaterFraction(amrex::Real x, amrex::Real y, amrex::Real dx, amrex::Real dy, Parm const& parm)
{
    const amrex::Real sn = std::sin(parm.theta_w);
    const amrex::Real cs = std::cos(parm.theta_w);
    const amrex::Real a  = (x + 0.5 * dx - parm.x0) * sn - (y - 0.5 * dy) * cs;
    const amrex::Real p  = dx * sn;
    const amrex::Real q  = dy * cs;

    auto F = [] (amrex::Real s) -> amrex::Real { return (s > 0.0) ? s * s : 0.0; };
    const amrex::Real area = (F(a) - F(a - p) - F(a - q) + F(a - p - q)) / (2.0 * sn * cs);
    return amrex::min(amrex::max(area / (dx * dy), 0.0), 1.0);
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
    sharp = true;
    amrex::ignore_unused(z, AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    // ---------------------------------------
    // Inclined interface: the fraction of the cell in the water, between eps/2 and 1 - eps/2
    // ---------------------------------------
    const amrex::Real eps = 1.0e-8;

    const amrex::Real Alpha2 = 0.5 * eps + (1.0 - eps) * WaterFraction(x, y, dx, dy, parm);
    const amrex::Real Alpha1 = 1.0 - Alpha2;

    // ---------------------------------------
    // Fluids at rest, and the post-shock air behind the shock
    // ---------------------------------------
    amrex::Real rho1 = parm.rho_air;
    amrex::Real rho2 = parm.rho_water;
    amrex::Real u    = 0.0;
    amrex::Real v    = 0.0;
    amrex::Real w    = 0.0;
    amrex::Real P    = parm.P0;

    if (x < parm.x_s) {
        rho1 = parm.rho_ps;
        u    = parm.u_ps;
        P    = parm.P_ps;
    }

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
// Boundary conditions
//
// With run.lo_bc = 1 in x, the ghost cells on the left are filled with the post-shock air,
// flowing in at the velocity behind the shock, so the shock keeps its speed. The other
// boundaries use the standard codes (slip wall below, outflow elsewhere).
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
    amrex::ignore_unused(time);
    const auto lo  = amrex::lbound(geom.domain);
    const bool xlo = (iv[0] < lo.x) && (bcr[bcomp].lo(0) == BCType::ext_dir);
    if (!xlo) { return; }

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

    // Pure post-shock air, with the same trace of water as the initial condition far from the interface
    const Real Alpha2 = 0.5e-8;
    Array<Real, NSTATE> Primitive, Conservative;
    State2Primitive(Primitive, (1.0 - Alpha2) * parm->rho_ps, Alpha2 * parm->rho_water,
       AMREX_D_DECL(parm->u_ps, 0.0, 0.0), parm->P_ps, 1.0 - Alpha2);
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
// N coarse steps. Each call appends one line to <plot directory>/waves.csv with
//   x_shock  the incident shock in the air along the top row, the largest x where the pressure
//            in the air crosses (P0 + P_ps)/2
//   y_R      the refraction point, the largest y of an interface cell (0.1 < alpha_1 < 0.9)
//            with P > (P0 + P_ps)/2, the foot of the incident shock or of the Mach stem
//   S_t      the transmitted wave in the water along the bottom wall, the largest x where the
//            pressure in the water crosses P0 + 0.05 (P_ps - P0)
// The rows are taken on the composite grid, i.e. each cell not covered by a finer level, and
// the crossings are interpolated linearly between a cell and its +x neighbor.
// =======================================

// Pressure of cell (i,j,k) of the conservative state U
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real CellPressure(UserOutput::State const& U, int i, int j, int k, Parm const& parm)
{
    amrex::Array<Real, NSTATE> Cons, Prim;
    for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(i,j,k,n); }
    Conservative2Primitive(Prim, Cons, parm);
    Real a1r1, a2r2, a1, P;
    AMREX_D_TERM(Real ux;, Real uy;, Real uz;)
    Primitive2State(a1r1, a2r2, AMREX_D_DECL(ux, uy, uz), P, a1, Prim);
    return P;
}

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    Parm const& hparm  = HostParm(o);
    const Real  P_half = 0.5 * (hparm.P0 + hparm.P_ps);
    const Real  P_foot = hparm.P0 + 0.05 * (hparm.P_ps - hparm.P0);
    const Real  big    = 1.0e30;

    // One ghost cell, for the +x neighbor in the crossings
    const auto v = Measure(o, {Max, Max, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<3>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real alpha1 = U(i,j,k,INDEX_VolumeFraction1);
        const Real P      = CellPressure(U, i,   j, k, parm);
        const Real Pe     = CellPressure(U, i+1, j, k, parm);

        // Where the pressure falls through the level toward +x, by linear interpolation
        auto front = [&] (Real level) -> Real {
            if (P <= level) { return -big; }
            return (Pe <= level) ? c.x[0] + c.dx[0] * (P - level) / (P - Pe) : c.x[0];
        };

        const bool top    = (j == c.domain.bigEnd(1));
        const bool bottom = (j == c.domain.smallEnd(1));
        const Real x_shock = (top && alpha1 > 0.5) ? front(P_half) : -big;
        const Real y_R     = (alpha1 > 0.1 && alpha1 < 0.9 && P > P_half) ? c.x[1] : -big;
        const Real S_t     = (bottom && alpha1 < 0.5) ? front(P_foot) : -big;
        return {x_shock, y_R, S_t};
    });

    WriteRow(o, "waves.csv", {{"x_shock", v[0]}, {"y_R", v[1]}, {"S_t", v[2]}});
}

#endif  // IC_H_
