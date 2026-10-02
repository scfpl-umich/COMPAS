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
// Decay of a shear wave u = U sin(k (x_a - a_lo)), v = w = 0, with x_a the coordinate along the
// axis prob.axis (1: y, or 2: z in 3D), a_lo its lower bound, k = 2 pi prob.mode/L_a and L_a the
// length of the box along that axis, in a uniform mixture of the phases (volume fractions
// prob.alpha, densities prob.rho) at the uniform pressure prob.P0. With the mixture viscosity
// mu = sum_k alpha_k mu_k and the density rho, the amplitude decays as exp(-nu k^2 t),
// nu = mu/rho, as long as the Mach number U/c is small.
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
    amrex::ignore_unused(x, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    const amrex::Real k   = 2.0*PI*parm.mode/parm.L_a;
    const amrex::Real x_a = (parm.axis == 2) ? z : y;
    const amrex::Real u   = parm.U*std::sin(k*(x_a - parm.a_lo));
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only
    const amrex::Real P = parm.P0;

    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha    [i_Phase] = parm.alpha_phase[i_Phase];
        Alpha_rho[i_Phase] = Alpha[i_Phase]*parm.rho_phase[i_Phase];
        Alpha_P  [i_Phase] = Alpha[i_Phase]*P;
    }
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    P);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/shear.csv with the integrals over the composite grid of
//   mass_k             partial densities alpha_k rho_k
//   mom_x, mom_y       momentum
//   energy             total energy rho E
//   kinetic            kinetic energy rho |u|^2/2
//   ie_mix             rho E - rho |u|^2/2
//   ie_k               phase internal energies alpha_k rho_k e_k
// and the amplitude of the wave, amp = 2/V int u sin(k (x_a - a_lo)) dV, with V the volume of
// the box. The box is periodic, so mass, momentum and energy are conserved, while the kinetic
// energy goes into the internal energies, in the ratio alpha_k mu_k without relaxation.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    constexpr int N = 2*NPHASE + 6;
    const auto v = Measure<N>(o, Sum, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<N>
    {
        const int  i = c.i, j = c.j, k = c.k;
        Real rho = 0.0;
        for (int n = 0; n < NPHASE; ++n) { rho += U(i,j,k,INDEX_Mass1+n); }
        AMREX_D_TERM(const Real ru = U(i,j,k,INDEX_MomentumX);,
                     const Real rv = U(i,j,k,INDEX_MomentumY);,
                     const Real rw = U(i,j,k,INDEX_MomentumZ);)
        const Real m2 = AMREX_D_TERM(ru*ru, + rv*rv, + rw*rw);
        const Real ke = 0.5*m2/rho;
        const Real rE = U(i,j,k,INDEX_Total_Energy);
        const Real kw = 2.0*PI*parm.mode/parm.L_a;
#if (AMREX_SPACEDIM > 1)
        const Real my = rv;
#else
        const Real my = 0.0;
#endif
        Values<N> r;
        for (int n = 0; n < NPHASE; ++n) {
            r[n]          = U(i,j,k,INDEX_Mass1+n)*c.dV;
            r[NPHASE+6+n] = U(i,j,k,INDEX_Energy1+n)*c.dV;
        }
        r[NPHASE]   = ru*c.dV;
        r[NPHASE+1] = my*c.dV;
        r[NPHASE+2] = rE*c.dV;
        r[NPHASE+3] = ke*c.dV;
        r[NPHASE+4] = (rE - ke)*c.dV;
        r[NPHASE+5] = 2.0/parm.Volume*ru/rho*std::sin(kw*(c.x[AMREX_SPACEDIM > 1 ? parm.axis : 0] - parm.a_lo))*c.dV;
        return r;
    });

    std::vector<std::pair<std::string, Real>> columns;
    for (int n = 0; n < NPHASE; ++n) { columns.push_back({"mass_" + std::to_string(n+1), v[n]}); }
    columns.insert(columns.end(), {{"mom_x", v[NPHASE]}, {"mom_y", v[NPHASE+1]}, {"energy", v[NPHASE+2]},
                                   {"kinetic", v[NPHASE+3]}, {"ie_mix", v[NPHASE+4]}, {"amp", v[NPHASE+5]}});
    for (int n = 0; n < NPHASE; ++n) { columns.push_back({"ie_" + std::to_string(n+1), v[NPHASE+6+n]}); }
    WriteRow(o, "shear.csv", columns);
}

#endif  // IC_H_
