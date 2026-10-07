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
// PhaseDiffusion_FaceViscosity_D: the viscosity of a face from the mixture viscosities of its
// two cells, their harmonic mean 2 mu_m mu_c/(mu_m + mu_c) (0 if both are 0). Across a material
// interface parallel to the face the traction is continuous and the normal velocity gradient
// jumps, so the two half cells act in series: with the harmonic mean a layered shear flow
// across a sharp interface that lies on a cell face is exact, while the arithmetic mean
// overestimates the face viscosity and gives a first-order interface error. An interface inside
// a cell (a mixed cell, with the mixture viscosity) keeps a first-order error. The mean is written as mu_m (2 mu_c/(mu_m + mu_c)),
// which is mu_m exactly where mu_m = mu_c (a uniform mixture). A negative cell value (a mixture
// coefficient formed from a volume fraction outside [0, 1]) counts as 0: the harmonic mean of
// signed values has a pole at mu_m = -mu_c and would give a negative, unbounded face value.
// With the clamp the face value lies between 0 and 2 min(mu_m, mu_c), and non-negative cell
// values give the same face value as without it.
//
// PhaseDiffusion_DissipationFraction_D: the volume fraction with which a phase takes its share
// alpha_k Phi_k of the viscous dissipation, alpha_k where the phase has mass (alpha_k > 0 and
// alpha_k rho_k > 0) and 0 elsewhere, so a massless trace (a round-off alpha_k, or a negative one)
// receives no energy.
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
void PhaseDiffusion_FaceViscosity_D(amrex::Real& mu_Face,
                                    amrex::Real const& mu_m,
                                    amrex::Real const& mu_c)
{
    amrex::Real const mu_mp = (mu_m < 0.0) ? 0.0 : mu_m;
    amrex::Real const mu_cp = (mu_c < 0.0) ? 0.0 : mu_c;
    amrex::Real const Sum = mu_mp + mu_cp;
    mu_Face = (Sum > 0.0) ? mu_mp*(2.0*mu_cp/Sum) : 0.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real PhaseDiffusion_DissipationFraction_D(amrex::Real const& Alpha,
                                                 amrex::Real const& Alpha_rho)
{
    return (Alpha > 0.0 && Alpha_rho > 0.0) ? Alpha : 0.0;
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
