// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_
#define FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// Face quantities of FluxDiffusion from the cell values of FDM_Conservative2FluxDiffusion_D:
//   - the face velocity and velocity gradient (FDM_FaceVelocityGradient_D, V[b][s] the
//     velocity component b at the stencil point s);
//   - the face viscosities, the harmonic mean of the mixture viscosities of the cells m and c
//     (PhaseDiffusion_FaceViscosity_D);
//   - per phase, the temperature difference across the face, the harmonic face weight and
//     the thermal diffusivity bound for the time step (where the phase conducts).
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_FluxDiffusion_D(amrex::Array<amrex::Real,NSTATE>& Flux,
                         amrex::Real& dmax,
                         amrex::Array<amrex::Array<amrex::Real,FaceStencil::N>,AMREX_SPACEDIM> const& V,
                         amrex::Real const& mu_m,  amrex::Real const& mu_c,
                         amrex::Real const& mub_m, amrex::Real const& mub_c,
                         amrex::Real const& rho_m, amrex::Real const& rho_c,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_m,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_c,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_m,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_c,
                         amrex::Array<amrex::Real,NPHASE> const& T_m,
                         amrex::Array<amrex::Real,NPHASE> const& T_c,
                         amrex::Array<bool,NPHASE> const& Present_m,
                         amrex::Array<bool,NPHASE> const& Present_c,
            AMREX_D_DECL(amrex::Real const& dx,
                         amrex::Real const& dy,
                         amrex::Real const& dz),
                         Parm const& EOS_Parameter,
                         int ID_Dimension)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Grad_Vel;
    FDM_FaceVelocityGradient_D(Vel,
                               Grad_Vel,
                               V,
                  AMREX_D_DECL(dx, dy, dz),
                               ID_Dimension);

    amrex::Real mu, mub;
    PhaseDiffusion_FaceViscosity_D(mu, mu_m, mu_c);
    PhaseDiffusion_FaceViscosity_D(mub, mub_m, mub_c);

    amrex::Array<amrex::Real,AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx, dy, dz)};
    amrex::Array<amrex::Real,NPHASE> Weight, Grad_T;
    amrex::Real Chi_T = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        PhaseDiffusion_FaceWeight_D(Weight[i_Phase], Alpha_m[i_Phase], Alpha_c[i_Phase], Present_m[i_Phase], Present_c[i_Phase]);
        Grad_T[i_Phase] = 0.0;
        if (Weight[i_Phase] > 0.0){
            FDM_Gradient_D(Grad_T[i_Phase], T_m[i_Phase], T_c[i_Phase], dX[ID_Dimension]);
            // The face weight over the smaller volume fraction (between 1 and 2) bounds the
            // explicit limit of both cells
            amrex::Real rho_k = std::min(Alpha_rho_m[i_Phase]/Alpha_m[i_Phase], Alpha_rho_c[i_Phase]/Alpha_c[i_Phase]);
            amrex::Real Ratio = Weight[i_Phase]/std::min(Alpha_m[i_Phase], Alpha_c[i_Phase]);
            Chi_T = std::max(Chi_T, Ratio*EOS_Parameter.kappa[i_Phase]/(rho_k*EOS_Parameter.cv[i_Phase]));
        }
    }

    amrex::Real rho = std::min(rho_m, rho_c);

    FluxDiffusion(Flux,
                  dmax,
                  Vel,
                  Grad_Vel,
                  Weight,
                  Grad_T,
                  mu,
                  mub,
                  rho,
                  Chi_T,
                  EOS_Parameter,
                  ID_Dimension);
}

#endif
#endif
