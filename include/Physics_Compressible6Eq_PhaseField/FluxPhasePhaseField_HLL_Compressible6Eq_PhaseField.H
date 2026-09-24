// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_HLL_Compressible6Eq_PhaseField_H_
#define FluxPhasePhaseField_HLL_Compressible6Eq_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_HLL(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& FluxPhase,
                             amrex::Real& J1_rho1,
                             amrex::Real& J2_rho2,
                             amrex::Real& J1_rho1_EK1,
                             amrex::Real& J2_rho2_EK2,
                             amrex::Real& cmax,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseL,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseR,
                             amrex::Real const& Grad_Q,
                             Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> FluxPhaseL;
    amrex::Real J1_rho1L;
    amrex::Real J2_rho2L;
    amrex::Real J1_rho1_EK1L;
    amrex::Real J2_rho2_EK2L;
    amrex::Real VPhase1L;
    amrex::Real VPhase2L;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseL,
                                          J1_rho1L, J2_rho2L,
                                          J1_rho1_EK1L, J2_rho2_EK2L,
                                          VPhase1L, VPhase2L,
                                          ConservativePhaseL,
                                          Grad_Q,
                                          EOS_Parameter);
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> FluxPhaseR;
    amrex::Real J1_rho1R;
    amrex::Real J2_rho2R;
    amrex::Real J1_rho1_EK1R;
    amrex::Real J2_rho2_EK2R;
    amrex::Real VPhase1R;
    amrex::Real VPhase2R;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseR,
                                          J1_rho1R, J2_rho2R,
                                          J1_rho1_EK1R, J2_rho2_EK2R,
                                          VPhase1R, VPhase2R,
                                          ConservativePhaseR,
                                          Grad_Q,
                                          EOS_Parameter);
    
    amrex::Real Sp = std::max({ 0.0, -VPhase1L, -VPhase2L, -VPhase1R, -VPhase2R }) + 1.0e-10;
    amrex::Real Sm = std::min({ 0.0, -VPhase1L, -VPhase2L, -VPhase1R, -VPhase2R }) - 1.0e-10;
    
    for (int i_State = 0; i_State < (3+AMREX_SPACEDIM)*2; i_State++){
        FluxPhase[i_State] = ( Sp*FluxPhaseL[i_State] - Sm*FluxPhaseR[i_State] - Sp*Sm*(ConservativePhaseR[i_State] - ConservativePhaseL[i_State]) )/(Sp - Sm);
    }
    
    cmax = std::max( -Sm, Sp );

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
    AMREX_D_TERM(amrex::Real uL = (Alpha1_rho1_u1L + Alpha2_rho2_u2L)/(Alpha1_rho1L + Alpha2_rho2L);,
                 amrex::Real vL = (Alpha1_rho1_v1L + Alpha2_rho2_v2L)/(Alpha1_rho1L + Alpha2_rho2L);,
                 amrex::Real wL = (Alpha1_rho1_w1L + Alpha2_rho2_w2L)/(Alpha1_rho1L + Alpha2_rho2L) );
    amrex::Real EKL = 0.5*(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL));

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
    AMREX_D_TERM(amrex::Real uR = (Alpha1_rho1_u1R + Alpha2_rho2_u2R)/(Alpha1_rho1R + Alpha2_rho2R);,
                 amrex::Real vR = (Alpha1_rho1_v1R + Alpha2_rho2_v2R)/(Alpha1_rho1R + Alpha2_rho2R);,
                 amrex::Real wR = (Alpha1_rho1_w1R + Alpha2_rho2_w2R)/(Alpha1_rho1R + Alpha2_rho2R) );
    amrex::Real EKR = 0.5*(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR));

    J1_rho1 = ( Sp*J1_rho1L - Sm*J1_rho1R - Sp*Sm*(Alpha1_rho1R - Alpha1_rho1L) )/(Sp - Sm);
    J2_rho2 = ( Sp*J2_rho2L - Sm*J2_rho2R - Sp*Sm*(Alpha2_rho2R - Alpha2_rho2L) )/(Sp - Sm);
    J1_rho1_EK1 = ( Sp*J1_rho1_EK1L - Sm*J1_rho1_EK1R - Sp*Sm*(Alpha1_rho1R*EKR - Alpha1_rho1L*EKL) )/(Sp - Sm);
    J2_rho2_EK2 = ( Sp*J2_rho2_EK2L - Sm*J2_rho2_EK2R - Sp*Sm*(Alpha2_rho2R*EKR - Alpha2_rho2L*EKL) )/(Sp - Sm);
    
}

#endif