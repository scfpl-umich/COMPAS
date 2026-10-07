// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxDiffusion_Compressible6Eq_IE_NPhase_H_
#define FluxDiffusion_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// Viscous and conductive face flux of the N-phase six-equation model (internal energy form),
// face normal e_d (d = ID_Dimension). The core solves U_t = -div(F - F_D) with F_D this flux.
//
// The mixture stress tau = sum_k alpha_k tau_k, with tau_k = mu_k (grad u + grad u^T)
// + (muB_k - 2/3 mu_k)(div u) I, uses the face viscosities mu and mub, the harmonic means of
// the mixture viscosities sum_k alpha_k mu_k of the two cells (PhaseDiffusion_FaceViscosity_D):
//   momentum                       tau . e_d
//   phase internal energy k        q_k = w_k kappa_k dT_k/dx_d
//   total energy rho E             u_f . tau . e_d + sum_k q_k
// Each phase conducts with its own temperature and the harmonic face volume fraction w_k
// (PhaseDiffusion_FaceWeight_D). The heat flux of rho E is the sum of the phase heat fluxes,
// so the mixture energy is conservative. The viscous dissipation alpha_k tau_k : grad u of the
// phase internal energies is a cell source in FVM_SurfaceIntegral_NC.
//
// dmax, for the time step vnn dx^2/dmax, is a diffusivity: the larger of the momentum
// diffusivity (4/3 mu + mub)/rho, with rho the smaller of the two cell densities, and
// Chi_T, the largest thermal diffusivity of the phases that conduct through the face.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxDiffusion(amrex::Array<amrex::Real,NSTATE>& Flux,
                   amrex::Real& dmax,
                   amrex::Array<amrex::Real,AMREX_SPACEDIM> const& Vel,
                   amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Grad_Vel,
                   amrex::Array<amrex::Real,NPHASE> const& Weight,
                   amrex::Array<amrex::Real,NPHASE> const& Grad_T,
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

    for (int iState = 0; iState < NSTATE; ++iState){
        Flux[iState] = 0.0;
    }

    amrex::Real qC = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        amrex::Real qC_Phase = Weight[i_Phase]*EOS_Parameter.kappa[i_Phase]*Grad_T[i_Phase];
        Flux[INDEX_Energy1 + i_Phase] = qC_Phase;
        qC += qC_Phase;
    }

    for (int j = 0; j < AMREX_SPACEDIM; ++j){
        Flux[INDEX_MomentumX + j] = tau[j];
    }
    Flux[INDEX_MomentumX + ID_Dimension] += tau_b;
    Flux[INDEX_Total_Energy] = Work + qC;

    dmax = std::max((4.0/3.0*mu + mub)/rho, Chi_T);
}

#endif
#endif
