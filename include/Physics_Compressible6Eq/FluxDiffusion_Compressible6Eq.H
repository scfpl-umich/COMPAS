// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxDiffusion_Compressible6Eq_H_
#define FluxDiffusion_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS && DIFFUSION == true)

// =======================================
// Viscous and conductive face flux of the two-phase six-equation model, face normal e_d
// (d = ID_Dimension). The core solves U_t = -div(F - F_D) with F_D this flux.
//
// The mixture stress tau = sum_k alpha_k tau_k, with tau_k = mu_k (grad u + grad u^T)
// + (muB_k - 2/3 mu_k)(div u) I, uses the face viscosities mu and mub of the arithmetic face
// volume fractions, as in the five-equation models:
//   momentum   tau . e_d
//   energy 1   u_f . tau . e_d + w_1 kappa_1 dT_1/dx_d
//   energy 2                     w_2 kappa_2 dT_2/dx_d
// The whole viscous energy flux u_f . tau . e_d goes to phase 1. It is also returned in
// Flux_Work, the non-conservative face component INDEX_NC_ViscousWork, from which
// FVM_SurfaceIntegral_NC moves the share of phase 2 (see FVM_SurfaceIntegralNC_Compressible6Eq.H).
// Each phase conducts with its own temperature and the harmonic face volume fraction w_k
// (PhaseDiffusion_FaceWeight_D).
//
// dmax, for the time step vnn dx^2/dmax, is a diffusivity: the larger of the momentum
// diffusivity (4/3 mu + mub)/rho, with rho the smaller of the two cell densities, and
// Chi_T, the largest thermal diffusivity of the phases that conduct through the face.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxDiffusion(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Flux,
                   amrex::Real& Flux_Work,
                   amrex::Real& dmax,
                   amrex::Array<amrex::Real,AMREX_SPACEDIM> const& Vel,
                   amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Grad_Vel,
                   amrex::Array<amrex::Real,2> const& Weight,
                   amrex::Array<amrex::Real,2> const& Grad_T,
                   amrex::Real const& mu,
                   amrex::Real const& mub,
                   amrex::Real const& rho,
                   amrex::Real const& Chi_T,
                   Parm const& EOS_Parameter,
                   int ID_Dimension)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> tau;
    amrex::Real Work = 0.0;
    for (int j = 0; j < AMREX_SPACEDIM; ++j){
        tau[j] = mu*(Grad_Vel[j][ID_Dimension] + Grad_Vel[ID_Dimension][j]);
        Work += Vel[j]*tau[j];
    }
    amrex::Real Div_Vel = AMREX_D_TERM( Grad_Vel[0][0],
                                      + Grad_Vel[1][1],
                                      + Grad_Vel[2][2] );
    amrex::Real tau_b = (mub - 2.0/3.0*mu)*Div_Vel;
    Work += Vel[ID_Dimension]*tau_b;

    amrex::Real qC1 = Weight[0]*EOS_Parameter.kappa_1*Grad_T[0];
    amrex::Real qC2 = Weight[1]*EOS_Parameter.kappa_2*Grad_T[1];

    for (int iState = 0; iState < 5+AMREX_SPACEDIM; ++iState){
        Flux[iState] = 0.0;
    }
    for (int j = 0; j < AMREX_SPACEDIM; ++j){
        Flux[INDEX_MomentumX+j] = tau[j];
    }
    Flux[INDEX_MomentumX+ID_Dimension] += tau_b;
    Flux[INDEX_Energy1] = Work + qC1;
    Flux[INDEX_Energy2] = qC2;

    Flux_Work = Work;

    dmax = std::max((4.0/3.0*mu + mub)/rho, Chi_T);
}

#endif
#endif
