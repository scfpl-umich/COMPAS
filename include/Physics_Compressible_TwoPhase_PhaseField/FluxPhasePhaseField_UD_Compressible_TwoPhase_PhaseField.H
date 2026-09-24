// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_UD_Compressible_TwoPhase_PhaseField_H_
#define FluxPhasePhaseField_UD_Compressible_TwoPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_UD(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& FluxPhase,
                            amrex::Real& cmax,
                            amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseL,
                            amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseR,
                            amrex::Real const& Grad_Q,
                            Parm const& EOS_Parameter)
{
    amrex::Real Grad_Q1 =  Grad_Q;
    amrex::Real Grad_Q2 = -Grad_Q1;
    
    amrex::Real Alpha1_rho1L;
    amrex::Real Alpha2_rho2L;
    amrex::Real Alpha1_rho1_u1L;
    amrex::Real Alpha2_rho2_u2L;
#if (AMREX_SPACEDIM > 1)
    amrex::Real Alpha1_rho1_v1L;
    amrex::Real Alpha2_rho2_v2L;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Alpha1_rho1_w1L;
    amrex::Real Alpha2_rho2_w2L;
#endif
    amrex::Real Alpha1_rho1_E1L;
    amrex::Real Alpha2_rho2_E2L;
    amrex::Real Alpha1L;
    amrex::Real Alpha2L;
    ConservativePhase2State(Alpha1_rho1L, Alpha2_rho2L,
                            Alpha1_rho1_u1L, Alpha2_rho2_u2L,
#if (AMREX_SPACEDIM > 1)
                            Alpha1_rho1_v1L, Alpha2_rho2_v2L,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha1_rho1_w1L, Alpha2_rho2_w2L,
#endif
                            Alpha1_rho1_E1L, Alpha2_rho2_E2L,
                            Alpha1L, Alpha2L,
                            ConservativePhaseL);
    
    amrex::Real Alpha1_rho1R;
    amrex::Real Alpha2_rho2R;
    amrex::Real Alpha1_rho1_u1R;
    amrex::Real Alpha2_rho2_u2R;
#if (AMREX_SPACEDIM > 1)
    amrex::Real Alpha1_rho1_v1R;
    amrex::Real Alpha2_rho2_v2R;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Alpha1_rho1_w1R;
    amrex::Real Alpha2_rho2_w2R;
#endif
    amrex::Real Alpha1_rho1_E1R;
    amrex::Real Alpha2_rho2_E2R;
    amrex::Real Alpha1R;
    amrex::Real Alpha2R;
    ConservativePhase2State(Alpha1_rho1R, Alpha2_rho2R,
                            Alpha1_rho1_u1R, Alpha2_rho2_u2R,
#if (AMREX_SPACEDIM > 1)
                            Alpha1_rho1_v1R, Alpha2_rho2_v2R,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha1_rho1_w1R, Alpha2_rho2_w2R,
#endif
                            Alpha1_rho1_E1R, Alpha2_rho2_E2R,
                            Alpha1R, Alpha2R,
                            ConservativePhaseR);
    
    //cmax
    cmax = 0.0;
    amrex::Real Alpha2_UGradQ1;
    f_Upwind(Alpha2_UGradQ1,
             Alpha2L,
             Alpha2R,
             Grad_Q1);
    cmax = std::max( cmax,std::abs(Alpha2_UGradQ1*Grad_Q1) );
    
    amrex::Real Alpha1_UGradQ2;
    f_Upwind(Alpha1_UGradQ2,
             Alpha1L,
             Alpha1R,
             Grad_Q2);
    cmax = std::max( cmax,std::abs(Alpha1_UGradQ2*Grad_Q2) );
    //Mass
    amrex::Real Alpha1_rho1_DGradQ1;
    f_Downwind(Alpha1_rho1_DGradQ1,
               Alpha1_rho1L,
               Alpha1_rho1R,
               Grad_Q1);
    amrex::Real J1_rho1 = Alpha1_rho1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_rho2_DGradQ2;
    f_Downwind(Alpha2_rho2_DGradQ2,
               Alpha2_rho2L,
               Alpha2_rho2R,
               Grad_Q2);
    amrex::Real J2_rho2 = Alpha2_rho2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
    //Momentum_x
    amrex::Real Alpha1_rho1_u1_DGradQ1;
    f_Downwind(Alpha1_rho1_u1_DGradQ1,
               Alpha1_rho1_u1L,
               Alpha1_rho1_u1R,
               Grad_Q1);
    amrex::Real J1_rho1_u1 = Alpha1_rho1_u1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_rho2_u2_DGradQ2;
    f_Downwind(Alpha2_rho2_u2_DGradQ2,
               Alpha2_rho2_u2L,
               Alpha2_rho2_u2R,
               Grad_Q2);
    amrex::Real J2_rho2_u2 = Alpha2_rho2_u2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
#if (AMREX_SPACEDIM > 1)
    //Momentum_y
    amrex::Real Alpha1_rho1_v1_DGradQ1;
    f_Downwind(Alpha1_rho1_v1_DGradQ1,
               Alpha1_rho1_v1L,
               Alpha1_rho1_v1R,
               Grad_Q1);
    amrex::Real J1_rho1_v1 = Alpha1_rho1_v1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_rho2_v2_DGradQ2;
    f_Downwind(Alpha2_rho2_v2_DGradQ2,
               Alpha2_rho2_v2L,
               Alpha2_rho2_v2R,
               Grad_Q2);
    amrex::Real J2_rho2_v2 = Alpha2_rho2_v2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
#endif
#if (AMREX_SPACEDIM > 2)
    //Momentum_z
    amrex::Real Alpha1_rho1_w1_DGradQ1;
    f_Downwind(Alpha1_rho1_w1_DGradQ1,
               Alpha1_rho1_w1L,
               Alpha1_rho1_w1R,
               Grad_Q1);
    amrex::Real J1_rho1_w1 = Alpha1_rho1_w1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_rho2_w2_DGradQ2;
    f_Downwind(Alpha2_rho2_w2_DGradQ2,
               Alpha2_rho2_w2L,
               Alpha2_rho2_w2R,
               Grad_Q2);
    amrex::Real J2_rho2_w2 = Alpha2_rho2_w2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
#endif
    //Energy
    amrex::Real Alpha1_rho1_E1_DGradQ1;
    f_Downwind(Alpha1_rho1_E1_DGradQ1,
               Alpha1_rho1_E1L,
               Alpha1_rho1_E1R,
               Grad_Q1);
    amrex::Real J1_rho1_E1 = Alpha1_rho1_E1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_rho2_E2_DGradQ2;
    f_Downwind(Alpha2_rho2_E2_DGradQ2,
               Alpha2_rho2_E2L,
               Alpha2_rho2_E2R,
               Grad_Q2);
    amrex::Real J2_rho2_E2 = Alpha2_rho2_E2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
    //Alpha1
    amrex::Real Alpha1_DGradQ1;
    f_Downwind(Alpha1_DGradQ1,
               Alpha1L,
               Alpha1R,
               Grad_Q1);
    amrex::Real J1 = Alpha1_DGradQ1*Alpha2_UGradQ1*Grad_Q1;
    
    amrex::Real Alpha2_DGradQ2;
    f_Downwind(Alpha2_DGradQ2,
               Alpha2L,
               Alpha2R,
               Grad_Q2);
    amrex::Real J2 = Alpha2_DGradQ2*Alpha1_UGradQ2*Grad_Q2;
    
    Component2VectorPhase(FluxPhase,
                          J1_rho1, J2_rho2,
                          J1_rho1_u1, J2_rho2_u2,
#if (AMREX_SPACEDIM > 1)
                          J1_rho1_v1, J2_rho2_v2,
#endif
#if (AMREX_SPACEDIM > 2)
                          J1_rho1_w1, J2_rho2_w2,
#endif
                          J1_rho1_E1, J2_rho2_E2,
                          J1, J2);
    
}

#endif