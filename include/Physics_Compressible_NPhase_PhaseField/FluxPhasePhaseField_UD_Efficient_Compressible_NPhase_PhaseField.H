// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_UD_Efficient_Compressible_NPhase_PhaseField_H_
#define FluxPhasePhaseField_UD_Efficient_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_UD_Efficient(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& FluxPhase,
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
    amrex::Real Sum_J = 0.0;
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
        
        amrex::Real AlphajL = 1.0 - AlphaiL;
        amrex::Real AlphajR = 1.0 - AlphaiR;
        
        //cmax
        amrex::Real Alphaj_UGradQi;
        f_Upwind(Alphaj_UGradQi,
                 AlphajL,
                 AlphajR,
                 Grad_Qi);
        amrex::Real cmaxi = std::abs(Alphaj_UGradQi*Grad_Qi);
        cmax = std::max( cmax,cmaxi );
        //Mass
        amrex::Real Alphai_rhoi_DGradQi;
        f_Downwind(Alphai_rhoi_DGradQi,
                   Alphai_rhoiL,
                   Alphai_rhoiR,
                   Grad_Qi);
        J_rho[i_Phase] = Alphai_rhoi_DGradQi*Alphaj_UGradQi*Grad_Qi;
        //Momentum_x
        amrex::Real Alphai_rhoi_ui_DGradQi;
        f_Downwind(Alphai_rhoi_ui_DGradQi,
                   Alphai_rhoi_uiL,
                   Alphai_rhoi_uiR,
                   Grad_Qi);
        J_rho_u[i_Phase] = Alphai_rhoi_ui_DGradQi*Alphaj_UGradQi*Grad_Qi;
#if (AMREX_SPACEDIM > 1)
        //Momentum_y
        amrex::Real Alphai_rhoi_vi_DGradQi;
        f_Downwind(Alphai_rhoi_vi_DGradQi,
                   Alphai_rhoi_viL,
                   Alphai_rhoi_viR,
                   Grad_Qi);
        J_rho_v[i_Phase] = Alphai_rhoi_vi_DGradQi*Alphaj_UGradQi*Grad_Qi;
#endif
#if (AMREX_SPACEDIM > 2)
        //Momentum_z
        amrex::Real Alphai_rhoi_wi_DGradQi;
        f_Downwind(Alphai_rhoi_wi_DGradQi,
                   Alphai_rhoi_wiL,
                   Alphai_rhoi_wiR,
                   Grad_Qi);
        J_rho_w[i_Phase] = Alphai_rhoi_wi_DGradQi*Alphaj_UGradQi*Grad_Qi;
#endif
        //Energy
        amrex::Real Alphai_rhoi_Ei_DGradQi;
        f_Downwind(Alphai_rhoi_Ei_DGradQi,
                   Alphai_rhoi_EiL,
                   Alphai_rhoi_EiR,
                   Grad_Qi);
        J_rho_E[i_Phase] = Alphai_rhoi_Ei_DGradQi*Alphaj_UGradQi*Grad_Qi;
        //Alpha
        amrex::Real Alphai_DGradQi;
        f_Downwind(Alphai_DGradQi,
                   AlphaiL,
                   AlphaiR,
                   Grad_Qi);
        J[i_Phase] = Alphai_DGradQi*Alphaj_UGradQi*Grad_Qi;
        //Sum_J
        Sum_J += J[i_Phase];
    }
    cmax += std::abs(Sum_J);
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        //Mass
        amrex::Real Alphai_rhoiL = Alpha_rhoL[i_Phase];
        amrex::Real Alphai_rhoiR = Alpha_rhoR[i_Phase];
        amrex::Real Alphai_rhoi_USumJ;
        f_Upwind(Alphai_rhoi_USumJ,
                 Alphai_rhoiL,
                 Alphai_rhoiR,
                 Sum_J);
        J_rho[i_Phase] += - Alphai_rhoi_USumJ*Sum_J;
        //Momentum_x
        amrex::Real Alphai_rhoi_uiL = Alpha_rho_uL[i_Phase];
        amrex::Real Alphai_rhoi_uiR = Alpha_rho_uR[i_Phase];
        amrex::Real Alphai_rhoi_ui_USumJ;
        f_Upwind(Alphai_rhoi_ui_USumJ,
                 Alphai_rhoi_uiL,
                 Alphai_rhoi_uiR,
                 Sum_J);
        J_rho_u[i_Phase] += - Alphai_rhoi_ui_USumJ*Sum_J;
#if (AMREX_SPACEDIM > 1)
        //Momentum_y
        amrex::Real Alphai_rhoi_viL = Alpha_rho_vL[i_Phase];
        amrex::Real Alphai_rhoi_viR = Alpha_rho_vR[i_Phase];
        amrex::Real Alphai_rhoi_vi_USumJ;
        f_Upwind(Alphai_rhoi_vi_USumJ,
                 Alphai_rhoi_viL,
                 Alphai_rhoi_viR,
                 Sum_J);
        J_rho_v[i_Phase] += - Alphai_rhoi_vi_USumJ*Sum_J;
#endif
#if (AMREX_SPACEDIM > 2)
        //Momentum_z
        amrex::Real Alphai_rhoi_wiL = Alpha_rho_wL[i_Phase];
        amrex::Real Alphai_rhoi_wiR = Alpha_rho_wR[i_Phase];
        amrex::Real Alphai_rhoi_wi_USumJ;
        f_Upwind(Alphai_rhoi_wi_USumJ,
                 Alphai_rhoi_wiL,
                 Alphai_rhoi_wiR,
                 Sum_J);
        J_rho_w[i_Phase] += - Alphai_rhoi_wi_USumJ*Sum_J;
#endif
        //Energy
        amrex::Real Alphai_rhoi_EiL = Alpha_rho_EL[i_Phase];
        amrex::Real Alphai_rhoi_EiR = Alpha_rho_ER[i_Phase];
        amrex::Real Alphai_rhoi_Ei_USumJ;
        f_Upwind(Alphai_rhoi_Ei_USumJ,
                 Alphai_rhoi_EiL,
                 Alphai_rhoi_EiR,
                 Sum_J);
        J_rho_E[i_Phase] += - Alphai_rhoi_Ei_USumJ*Sum_J;
        //Alpha
        amrex::Real AlphaiL = AlphaL[i_Phase];
        amrex::Real AlphaiR = AlphaR[i_Phase];
        amrex::Real Alphai_USumJ;
        f_Upwind(Alphai_USumJ,
                 AlphaiL,
                 AlphaiR,
                 Sum_J);
        J[i_Phase] += - Alphai_USumJ*Sum_J;
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