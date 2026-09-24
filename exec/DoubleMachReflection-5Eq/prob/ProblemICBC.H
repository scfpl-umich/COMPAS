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
// Double Mach reflection (Woodward and Colella, J. Comput. Phys. 54, 1984), with the
// boundary treatment discussed by Kemm (Comput. Fluids 132, 2016).
// A Mach 10 shock in air runs up a 30 degree wedge. The frame is aligned with the wedge,
// so its surface is the bottom boundary for x >= x0 = 1/6, and the shock is inclined at
// 60 degrees to it. The shock passes through (x0, 0) at t = 0 and moves along its normal
// (cos theta_w, -sin theta_w) at the speed Us, so at height y and time t it is at
//   x_s(y, t) = x0 + y tan(theta_w) + Us t / cos(theta_w).
// Single fluid, gamma = 1.4, in the five-equation model with alpha_1 = 1.
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real shock_position(amrex::Real y, amrex::Real t, Parm const& parm)
{
    return parm.x0 + y * parm.tan_w + parm.Us * t / parm.cos_w;
}


// Conservative state behind (post = true) or ahead of the shock
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void shock_state(amrex::Array<amrex::Real, NSTATE>& Conservative, bool post, Parm const& parm)
{
    const amrex::Real rho = post ? parm.rhoS : parm.rho0;
    const amrex::Real u   = post ? parm.uS   : 0.0;
    const amrex::Real v   = post ? parm.vS   : 0.0;
    const amrex::Real w   = 0.0;
    const amrex::Real P   = post ? parm.PS   : parm.P0;

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
// Initial condition
//
// The state behind the shock to the left of x_s(y, 0), gas at rest to the right. The cell
// centers are sampled, so the oblique shock starts as a staircase on the grid. Kemm (2016)
// shows that smearing it does not remove the start-up disturbance that follows the shock.
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
    amrex::ignore_unused(AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    shock_state(Conservative, x < shock_position(y, 0.0, parm), parm);
}


// =======================================
// Optional: boundary conditions
//
// With run.lo_bc = 1 1 and run.hi_bc = 2 1, this function fills the ghost cells of the
// left, bottom and top faces (the right face is an outflow):
//   left    the state behind the shock
//   bottom  the state behind the shock ahead of the wedge (x < x0), and a slip wall on the
//           wedge (x >= x0): the mirror image of the interior with rho v negated
//   top     the exact moving shock, i.e. the state behind the shock for x < x_s(y, time)
//           and the gas at rest beyond, at the center of each ghost cell
// In the corner below the outflow face, the mirror source is clamped into the domain.
// The top boundary follows the undisturbed shock, which is exact until the reflected shock
// reaches y = 1. Kemm (2016) shows that the under-resolved shock in these ghost cells sends a
// weak disturbance into the domain, and moves it away from the triple points with the taller
// domain [0, 4] x [0, 2] (geometry.prob_hi = 4 2, amr.n_cell = 240 120).
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
    const auto lo  = amrex::lbound(geom.domain);
    const auto hi  = amrex::ubound(geom.domain);
    const bool xlo = (iv[0] < lo.x) && (bcr[bcomp].lo(0) == BCType::ext_dir);
    const bool ylo = (iv[1] < lo.y) && (bcr[bcomp].lo(1) == BCType::ext_dir);
    const bool yhi = (iv[1] > hi.y) && (bcr[bcomp].hi(1) == BCType::ext_dir);
    if (!xlo && !ylo && !yhi) { return; }

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

    // Slip wall on the wedge: mirror the interior across y = 0 and reverse rho v. The source
    // column is clamped so the corner cells below the outflow face read from the domain.
    if (ylo && x >= parm->x0) {
        IntVect im = iv;
        im[0] = amrex::max(lo.x, amrex::min(iv[0], hi.x));
        im[1] = 2 * lo.y - iv[1] - 1;
        for (int n = dcomp; n < dcomp + numcomp; ++n) {
            const int s = orig_comp + n - dcomp;
            q(iv, n) = (s == INDEX_MomentumY) ? -q(im, n) : q(im, n);
        }
        return;
    }

    // Left face, the bottom ahead of the wedge, and the top face: the exact shock solution
    const bool post = ylo || (xlo && !yhi) || (x < shock_position(y, time, *parm));
    Array<Real, NSTATE> Conservative;
    shock_state(Conservative, post, *parm);
    for (int n = dcomp; n < dcomp + numcomp; ++n) {
        const int s = orig_comp + n - dcomp;
        if (s < NSTATE) { q(iv, n) = Conservative[s]; }
    }
}


// =======================================
// Optional: in-situ diagnostics
//
// With run.user_output_int = N > 0 in the inputs file, COMPAS calls this function every
// N coarse steps. Each call appends one line to <plot directory>/stem.csv with
//   x_stem    foot of the Mach stem on the wall, the furthest point in x along the first
//             row of cells on the finest level where the pressure crosses halfway between
//             the pressures ahead of and behind the incident shock
//   x_inc     where the undisturbed incident shock would meet the wall, x_s(0, t)
// In a self-similar reflection x_stem - x0 grows linearly in t, faster than x_inc - x0.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    Parm const& hparm = HostParm(o);
    const Real  P_thr = 0.5 * (hparm.P0 + hparm.PS);

    // Finest level only, with one ghost cell for the crossing with the +x neighbor
    const auto v = Measure(o, {Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<1>
    {
        if (c.j != c.domain.smallEnd(1)) { return {-1.0e30}; }
        amrex::Array<Real, NSTATE> Uc, Ue;
        for (int n = 0; n < NSTATE; ++n) { Uc[n] = U(c.i,c.j,c.k,n); Ue[n] = U(c.i+1,c.j,c.k,n); }
        Real Pc, Pe;
        Conservative2Pressure(Pc, Uc, parm);
        Conservative2Pressure(Pe, Ue, parm);
        const bool cross = (Pc - P_thr) * (Pe - P_thr) <= 0.0 && Pc != Pe;
        return {cross ? c.x[0] + c.dx[0] * (P_thr - Pc) / (Pe - Pc) : -1.0e30};
    }, FinestLevel);
    WriteRow(o, "stem.csv", {{"x_stem", v[0]}, {"x_inc", shock_position(0.0, Time(o), hparm)}});
}

#endif  // IC_H_
