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
// Compressible Kelvin-Helmholtz instability of the shear layer u = U tanh(y/delta) at uniform
// density and pressure (Blumen, J. Fluid Mech. 40, 1970). The pressure P0 = rho0 U^2/(gamma M^2)
// sets Blumen's Mach number M = U/a. Both materials are the same gas, so the volume fraction
// alpha_1 = (1 + tanh(y/(w delta)))/2, with w = prob.tracer_width (1 by default), is a passive
// tracer of the layer; a small w gives a sharp tracer interface at y = 0. The layer is seeded with one
// divergence-free mode of wavenumber kx = 2 pi/Lx, with stream function
//   psi = (eps U/kx) cos(kx x) sech(y/delta),  u' = dpsi/dy,  v' = -dpsi/dx = eps U sin(kx x) sech(y/delta)
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
    amrex::ignore_unused(z, AMREX_D_DECL(dx, dy, dz), AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    // ---------------------------------------
    // Passive tracer: material 1 above the layer, material 2 below
    // ---------------------------------------
    const amrex::Real eta    = y / parm.delta;
    const amrex::Real eta_a  = y / (parm.tracer_width * parm.delta);    // equals eta for the default width 1
    const amrex::Real eps_a  = 1.0e-8;
    const amrex::Real Alpha1 = (parm.tracer == 1) ? 0.5 * (1.0 + std::tanh(eta_a) * (1.0 - eps_a)) : 1.0 - 0.5 * eps_a;
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    // ---------------------------------------
    // Shear layer and the seeded mode
    // ---------------------------------------
    const amrex::Real sech = 1.0 / std::cosh(eta);
    const amrex::Real u = parm.U * std::tanh(eta)
                        - parm.eps * parm.U / (parm.kx * parm.delta) * std::cos(parm.kx * x) * sech * std::tanh(eta);
    const amrex::Real v = parm.eps * parm.U * std::sin(parm.kx * x) * sech;
    const amrex::Real w = 0.0;
    const amrex::Real P = parm.P0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1 * parm.rho0, Alpha2 * parm.rho0,
       AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Spatio-temporal refinement tagging
//
// Tags the band |y| < prob.y_ref[lev] on level lev. The band is the same at every x, so the
// grids do not break the symmetry of the base flow in x and do not seed other modes.
// =======================================

#define SPATIO_TEMPORAL_TAG_USER_FUNC_
#define SPACE_TIME_REF
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
spatio_temporal_tag (int & tagged,
                     int lev,
                     amrex::Real x,
                     amrex::Real y,
                     amrex::Real z,
                     amrex::Real t,
                     Parm const& parm)
{
    amrex::ignore_unused(x, z, t);
    if (lev < 4 && std::abs(y) < parm.y_ref[lev]){
        tagged = TagBox::SET;
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends one
// line to <plot directory>/mode.csv with
//   S, C     int v sin(kx x) dA and int v cos(kx x) dA over the domain
//   A        amplitude of the seeded mode, 2 sqrt(S^2 + C^2)/Lx, which is eps U int sech dy at t = 0
//   Ey       y kinetic energy, int rho v^2/2 dA
//   Ey_max   largest rho v^2/2
//   theta    momentum thickness of the layer, int (U^2 - u^2)/(4 U^2) dA / Lx (delta/2 for the tanh profile)
// The integrals are taken over the composite grid, i.e. each cell not covered by a finer level.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    Parm const& hparm = HostParm(o);
    const Real  Lx    = 2.0 * PI / hparm.kx;

    const auto m = Measure(o, {Sum, Sum, Sum, Max, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        const Real rho = U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2);
        const Real u   = U(i,j,k,INDEX_MomentumX) / rho;
        const Real v   = U(i,j,k,INDEX_MomentumY) / rho;
        const Real ey  = 0.5 * rho * v * v;
        const Real th  = 0.25 * (1.0 - u * u / (parm.U * parm.U));
        return {v * std::sin(parm.kx * c.x[0]) * c.dV, v * std::cos(parm.kx * c.x[0]) * c.dV, ey * c.dV, ey, th * c.dV};
    });

    WriteRow(o, "mode.csv", {{"S", m[0]}, {"C", m[1]}, {"A", 2.0 * std::sqrt(m[0]*m[0] + m[1]*m[1]) / Lx},
                             {"Ey", m[2]}, {"Ey_max", m[3]}, {"theta", m[4] / Lx}});
}

#endif  // IC_H_
