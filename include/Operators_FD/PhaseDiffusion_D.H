// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseDiffusion_D_H_
#define PhaseDiffusion_D_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

// =======================================
// Pieces of the viscous and conductive fluxes of the six-equation models (DIFFUSION = true)
// that do not depend on the model.
//
// Phase k is present in a cell where alpha_k >= PhaseDiffusion::AlphaMin and alpha_k rho_k > 0.
// Its temperature is defined there only.
//
// PhaseDiffusion_FaceWeight_D: the volume fraction that carries the heat flux of phase k
// through a face, the harmonic mean w = 2 a_m a_c/(a_m + a_c) of the two cells, and 0 unless
// the phase is present in both. A phase conducts only through itself, so an absent phase
// receives no heat, and since w <= 2 min(a_m, a_c) the explicit stability limit stays bounded
// by the phase's own heat capacity.
//
// ViscousDissipation_D: the dissipation tau : grad u of a phase with the stress
//   tau = mu (grad u + grad u^T) + (muB - 2/3 mu)(div u) I,
//   Phi = mu (d_i u_j + d_j u_i) d_j u_i + (muB - 2/3 mu)(div u)^2
//       = 2 mu |S_dev|^2 + muB (tr S)^2 >= 0     (S the symmetric part of grad u),
// with Grad_Vel[a][b] = d u_b / d x_a.
// =======================================

namespace PhaseDiffusion
{
constexpr amrex::Real AlphaMin = 1.0e-10;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool PhaseDiffusion_Present(amrex::Real const& Alpha,
                            amrex::Real const& Alpha_rho)
{
    return (Alpha >= PhaseDiffusion::AlphaMin) && (Alpha_rho > 0.0);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void PhaseDiffusion_FaceWeight_D(amrex::Real& Weight,
                                 amrex::Real const& Alpha_m,
                                 amrex::Real const& Alpha_c,
                                 bool Present_m,
                                 bool Present_c)
{
    Weight = (Present_m && Present_c) ? 2.0*Alpha_m*Alpha_c/(Alpha_m + Alpha_c) : 0.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void ViscousDissipation_D(amrex::Real& Phi,
                          amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Grad_Vel,
                          amrex::Real const& mu,
                          amrex::Real const& mub)
{
    amrex::Real Shear = 0.0;
    for (int a = 0; a < AMREX_SPACEDIM; ++a){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Shear += (Grad_Vel[a][b] + Grad_Vel[b][a])*Grad_Vel[b][a];
        }
    }
    amrex::Real Div_Vel = AMREX_D_TERM( Grad_Vel[0][0],
                                      + Grad_Vel[1][1],
                                      + Grad_Vel[2][2] );
    Phi = mu*Shear + (mub - 2.0/3.0*mu)*Div_Vel*Div_Vel;
}

#endif
