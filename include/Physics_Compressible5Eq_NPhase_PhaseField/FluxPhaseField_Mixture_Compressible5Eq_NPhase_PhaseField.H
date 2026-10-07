// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Mixture_Compressible5Eq_NPhase_PhaseField_H_
#define FluxPhaseField_Mixture_Compressible5Eq_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhaseField_Mixture(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& ConservativeL,
                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& ConservativeR,
                            amrex::Array<amrex::Real,NPHASE> const& J_rho,
                            amrex::Array<amrex::Real,NPHASE> const& J,
                            Parm const& EOS_Parameter)
{
    //Primitive
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;
    amrex::Array<amrex::Real,NPHASE> AlphaL;
    Primitive2State(Alpha_rhoL,
       AMREX_D_DECL(uL, vL, wL),
                    pL,
                    AlphaL,
                    PrimitiveL);
    
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;
    amrex::Array<amrex::Real,NPHASE> AlphaR;
    Primitive2State(Alpha_rhoR,
       AMREX_D_DECL(uR, vR, wR),
                    pR,
                    AlphaR,
                    PrimitiveR);
    
    //Momentum
    amrex::Real M = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        M += J_rho[i_Phase];
    }
    amrex::Real u_DM;
    f_Downwind(u_DM,
               uL,
               uR,
               M);
    amrex::Real M_u = M*u_DM;
#if (AMREX_SPACEDIM > 1)
    amrex::Real v_DM;
    f_Downwind(v_DM,
               vL,
               vR,
               M);
    amrex::Real M_v = M*v_DM;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real w_DM;
    f_Downwind(w_DM,
               wL,
               wR,
               M);
    amrex::Real M_w = M*w_DM;
#endif
    //Energy
    amrex::Real eK_DM = 0.5*( AMREX_D_TERM(u_DM*u_DM,+v_DM*v_DM,+w_DM*w_DM) );
    amrex::Real M_eK = M*eK_DM;
    
    amrex::Real J_rho_e = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real Ji = J[i_Phase];
        amrex::Real p_DJi;
        f_Downwind(p_DJi,
                   pL,
                   pR,
                   Ji);
        // NASG: the A and B terms carry the flux of phi_k = alpha_k - b_k alpha_k rho_k, J_k - b_k J_rho_k
        amrex::Real Ji_Phi = (EOS_Parameter.C[i_Phase] > 0.0) ? Ji - EOS_Parameter.C[i_Phase]*J_rho[i_Phase] : Ji;
        J_rho_e += EOS_Parameter.A[i_Phase]*Ji_Phi*p_DJi + EOS_Parameter.B[i_Phase]*Ji_Phi + EOS_Parameter.D[i_Phase]*J_rho[i_Phase];
    }
    
    //Flux_PhaseField
    Component2Vector(Flux,
                     J_rho,
        AMREX_D_DECL(M_u, M_v, M_w),
                     M_eK + J_rho_e,
                     J);
    
}

#endif