// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_UD_Base_Compressible_NPhase_PhaseField_H_
#define FluxPhasePhaseField_UD_Base_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_UD_Base(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& FluxPhase,
                                 amrex::Real& cmax,
                                 amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhaseL,
                                 amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhaseR,
                                 amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                                 Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_uL;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_vL;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_wL;
#endif
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_EL;
    amrex::Array<amrex::Real,NPHASE> AlphaL;
    ConservativePhase2State(Alpha_rhoL,
                            Alpha_rho_uL,
#if (AMREX_SPACEDIM > 1)
                            Alpha_rho_vL,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha_rho_wL,
#endif
                            Alpha_rho_EL,
                            AlphaL,
                            ConservativePhaseL);

    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_uR;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_vR;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_wR;
#endif
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_ER;
    amrex::Array<amrex::Real,NPHASE> AlphaR;
    ConservativePhase2State(Alpha_rhoR,
                            Alpha_rho_uR,
#if (AMREX_SPACEDIM > 1)
                            Alpha_rho_vR,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha_rho_wR,
#endif
                            Alpha_rho_ER,
                            AlphaR,
                            ConservativePhaseR);

    cmax = 0.0;
    amrex::Array<amrex::Real,NPHASE> J_rho;
    amrex::Array<amrex::Real,NPHASE> J_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> J_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> J_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> J_rho_E;
    amrex::Array<amrex::Real,NPHASE> J;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real Grad_Qi = Grad_Q[i_Phase];
        
        amrex::Real Alphai_rhoiL = Alpha_rhoL[i_Phase];
        amrex::Real Alphai_rhoiR = Alpha_rhoR[i_Phase];
        
        amrex::Real Alphai_rhoi_uiL = Alpha_rho_uL[i_Phase];
        amrex::Real Alphai_rhoi_uiR = Alpha_rho_uR[i_Phase];
#if (AMREX_SPACEDIM > 1)
        amrex::Real Alphai_rhoi_viL = Alpha_rho_vL[i_Phase];
        amrex::Real Alphai_rhoi_viR = Alpha_rho_vR[i_Phase];
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real Alphai_rhoi_wiL = Alpha_rho_wL[i_Phase];
        amrex::Real Alphai_rhoi_wiR = Alpha_rho_wR[i_Phase];
#endif
        amrex::Real Alphai_rhoi_EiL = Alpha_rho_EL[i_Phase];
        amrex::Real Alphai_rhoi_EiR = Alpha_rho_ER[i_Phase];
        
        amrex::Real AlphaiL = AlphaL[i_Phase];
        amrex::Real AlphaiR = AlphaR[i_Phase];
        
        amrex::Real cmaxi = 0.0;
        amrex::Real Ji_rhoi = 0.0;
        amrex::Real Ji_rhoi_ui = 0.0;
#if (AMREX_SPACEDIM > 1)
        amrex::Real Ji_rhoi_vi = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real Ji_rhoi_wi = 0.0;
#endif
        amrex::Real Ji_rhoi_Ei = 0.0;
        amrex::Real Ji = 0.0;
        for (int j_Phase = 0; j_Phase < NPHASE; j_Phase++){
            amrex::Real Grad_Qj = Grad_Q[j_Phase];
            
            amrex::Real AlphajL = AlphaL[j_Phase];
            amrex::Real AlphajR = AlphaR[j_Phase];
            
            amrex::Real Grad_Qij = Grad_Qi - Grad_Qj;
            
            //cmax
            amrex::Real Alphaj_UGradQij;
            f_Upwind(Alphaj_UGradQij,
                     AlphajL,
                     AlphajR,
                     Grad_Qij);
            cmaxi += std::abs(Alphaj_UGradQij*Grad_Qij);
            //Mass
            amrex::Real Alphai_rhoi_DGradQij;
            f_Downwind(Alphai_rhoi_DGradQij,
                       Alphai_rhoiL,
                       Alphai_rhoiR,
                       Grad_Qij);
            Ji_rhoi += Alphai_rhoi_DGradQij*Alphaj_UGradQij*Grad_Qij;
            //Momentum_x
            amrex::Real Alphai_rhoi_ui_DGradQij;
            f_Downwind(Alphai_rhoi_ui_DGradQij,
                       Alphai_rhoi_uiL,
                       Alphai_rhoi_uiR,
                       Grad_Qij);
            Ji_rhoi_ui += Alphai_rhoi_ui_DGradQij*Alphaj_UGradQij*Grad_Qij;
#if (AMREX_SPACEDIM > 1)
            //Momentum_y
            amrex::Real Alphai_rhoi_vi_DGradQij;
            f_Downwind(Alphai_rhoi_vi_DGradQij,
                       Alphai_rhoi_viL,
                       Alphai_rhoi_viR,
                       Grad_Qij);
            Ji_rhoi_vi += Alphai_rhoi_vi_DGradQij*Alphaj_UGradQij*Grad_Qij;
#endif
#if (AMREX_SPACEDIM > 2)
            //Momentum_z
            amrex::Real Alphai_rhoi_wi_DGradQij;
            f_Downwind(Alphai_rhoi_wi_DGradQij,
                       Alphai_rhoi_wiL,
                       Alphai_rhoi_wiR,
                       Grad_Qij);
            Ji_rhoi_wi += Alphai_rhoi_wi_DGradQij*Alphaj_UGradQij*Grad_Qij;
#endif
            //Energy
            amrex::Real Alphai_rhoi_Ei_DGradQij;
            f_Downwind(Alphai_rhoi_Ei_DGradQij,
                       Alphai_rhoi_EiL,
                       Alphai_rhoi_EiR,
                       Grad_Qij);
            Ji_rhoi_Ei += Alphai_rhoi_Ei_DGradQij*Alphaj_UGradQij*Grad_Qij;
            //Alpha
            amrex::Real Alphai_DGradQij;
            f_Downwind(Alphai_DGradQij,
                       AlphaiL,
                       AlphaiR,
                       Grad_Qij);
            Ji += Alphai_DGradQij*Alphaj_UGradQij*Grad_Qij;
        }
        cmax = std::max( cmax,cmaxi );
        J_rho[i_Phase] = Ji_rhoi;
        J_rho_u[i_Phase] = Ji_rhoi_ui;
#if (AMREX_SPACEDIM > 1)
        J_rho_v[i_Phase] = Ji_rhoi_vi;
#endif
#if (AMREX_SPACEDIM > 2)
        J_rho_w[i_Phase] = Ji_rhoi_wi;
#endif
        J_rho_E[i_Phase] = Ji_rhoi_Ei;
        J[i_Phase] = Ji;
    }

    Component2VectorPhase(FluxPhase,
                          J_rho,
                          J_rho_u,
#if (AMREX_SPACEDIM > 1)
                          J_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                          J_rho_w,
#endif
                          J_rho_E,
                          J);
    
}

#endif