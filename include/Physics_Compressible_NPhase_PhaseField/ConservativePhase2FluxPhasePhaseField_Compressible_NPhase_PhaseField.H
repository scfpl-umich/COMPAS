// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativePhase2FluxPhasePhaseField_Compressible_NPhase_PhaseField_H_
#define ConservativePhase2FluxPhasePhaseField_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativePhase2FluxPhasePhaseField(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& FluxPhase,
                                           amrex::Array<amrex::Real,NPHASE>& VPhase,
                                           amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhase,
                                           amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                                           Parm const& EOS_Parameter,
                                           int ID_Formulation)
{
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_E;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    ConservativePhase2State(Alpha_rho,
                            Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                            Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha_rho_w,
#endif
                            Alpha_rho_E,
                            Alpha,
                            ConservativePhase);

    //Base formulation
    if (ID_Formulation == 0){
        VelocityPhaseField_Base(VPhase,
                                Alpha,
                                Grad_Q);
    }
    //Efficient formulation
    else if (ID_Formulation == 1){
        VelocityPhaseField_Efficient(VPhase,
                                     Alpha,
                                     Grad_Q);
    }

    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho_E;
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Flux_Alpha_rho[i_Phase] = Alpha_rho[i_Phase]*VPhase[i_Phase];
        Flux_Alpha_rho_u[i_Phase] = Alpha_rho_u[i_Phase]*VPhase[i_Phase];
#if (AMREX_SPACEDIM > 1)
        Flux_Alpha_rho_v[i_Phase] = Alpha_rho_v[i_Phase]*VPhase[i_Phase];
#endif
#if (AMREX_SPACEDIM > 2)
        Flux_Alpha_rho_w[i_Phase] = Alpha_rho_w[i_Phase]*VPhase[i_Phase];
#endif
        Flux_Alpha_rho_E[i_Phase] = Alpha_rho_E[i_Phase]*VPhase[i_Phase];
        Flux_Alpha[i_Phase] = Alpha[i_Phase]*VPhase[i_Phase];
    }
    Component2VectorPhase(FluxPhase,
                          Flux_Alpha_rho,
                          Flux_Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                          Flux_Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                          Flux_Alpha_rho_w,
#endif
                          Flux_Alpha_rho_E,
                          Flux_Alpha);
    
}

#endif