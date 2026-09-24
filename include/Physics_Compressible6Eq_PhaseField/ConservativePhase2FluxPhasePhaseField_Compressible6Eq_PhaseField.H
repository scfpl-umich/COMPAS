// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativePhase2FluxPhasePhaseField_Compressible6Eq_PhaseField_H_
#define ConservativePhase2FluxPhasePhaseField_Compressible6Eq_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativePhase2FluxPhasePhaseField(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& FluxPhase,
                                           amrex::Real& J1_rho1,
                                           amrex::Real& J2_rho2,
                                           amrex::Real& J1_rho1_EK1,
                                           amrex::Real& J2_rho2_EK2,
                                           amrex::Real& VPhase1,
                                           amrex::Real& VPhase2,
                                           amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhase,
                                           amrex::Real const& Grad_Q,
                                           Parm const& EOS_Parameter)
{
    amrex::Real Grad_Q1 =  Grad_Q;
    amrex::Real Grad_Q2 = -Grad_Q1;
    
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    amrex::Real Alpha1_rho1_u1;
    amrex::Real Alpha2_rho2_u2;
#if (AMREX_SPACEDIM > 1)
    amrex::Real Alpha1_rho1_v1;
    amrex::Real Alpha2_rho2_v2;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Alpha1_rho1_w1;
    amrex::Real Alpha2_rho2_w2;
#endif
    amrex::Real Alpha1_rho1_E1;
    amrex::Real Alpha2_rho2_E2;
    amrex::Real Alpha1;
    amrex::Real Alpha2;
    ConservativePhase2State(Alpha1_rho1, Alpha2_rho2,
                            Alpha1_rho1_u1, Alpha2_rho2_u2,
#if (AMREX_SPACEDIM > 1)
                            Alpha1_rho1_v1, Alpha2_rho2_v2,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha1_rho1_w1, Alpha2_rho2_w2,
#endif
                            Alpha1_rho1_E1, Alpha2_rho2_E2,
                            Alpha1, Alpha2,
                            ConservativePhase);

    VPhase1 = Alpha2*Grad_Q1;
    VPhase2 = Alpha1*Grad_Q2;
    
    Component2VectorPhase(FluxPhase,
                          Alpha1_rho1*VPhase1, Alpha2_rho2*VPhase2,
                          Alpha1_rho1_u1*VPhase1, Alpha2_rho2_u2*VPhase2,
#if (AMREX_SPACEDIM > 1)
                          Alpha1_rho1_v1*VPhase1, Alpha2_rho2_v2*VPhase2,
#endif
#if (AMREX_SPACEDIM > 2)
                          Alpha1_rho1_w1*VPhase1, Alpha2_rho2_w2*VPhase2,
#endif
                          Alpha1_rho1_E1*VPhase1, Alpha2_rho2_E2*VPhase2,
                          Alpha1*VPhase1, Alpha2*VPhase2);

    AMREX_D_TERM(amrex::Real u = (Alpha1_rho1_u1 + Alpha2_rho2_u2)/(Alpha1_rho1 + Alpha2_rho2);,
                 amrex::Real v = (Alpha1_rho1_v1 + Alpha2_rho2_v2)/(Alpha1_rho1 + Alpha2_rho2);,
                 amrex::Real w = (Alpha1_rho1_w1 + Alpha2_rho2_w2)/(Alpha1_rho1 + Alpha2_rho2) );
    amrex::Real EK = 0.5*(AMREX_D_TERM(u*u, + v*v, + w*w));
    J1_rho1 = Alpha1_rho1*VPhase1;
    J2_rho2 = Alpha2_rho2*VPhase2;
    J1_rho1_EK1 = Alpha1_rho1*EK*VPhase1;
    J2_rho2_EK2 = Alpha2_rho2*EK*VPhase2;
    
}

#endif