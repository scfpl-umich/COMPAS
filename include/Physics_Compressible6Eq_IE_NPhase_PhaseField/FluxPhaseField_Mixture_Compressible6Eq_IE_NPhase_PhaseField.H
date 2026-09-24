// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Mixture_Compressible6Eq_IE_NPhase_PhaseField_H_
#define FluxPhaseField_Mixture_Compressible6Eq_IE_NPhase_PhaseField_H_

#if (PHYSICS == SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhaseField_Mixture(amrex::Array<amrex::Real, NSTATE>& Flux,
                            amrex::Array<amrex::Real, NSTATE> const& ConservativeL,
                            amrex::Array<amrex::Real, NSTATE> const& ConservativeR,
                            amrex::Array<amrex::Real,NPHASE> const& J_rho,
                            amrex::Array<amrex::Real,NPHASE> const& J,
                            Parm const& EOS_Parameter)
{
    //Primitive
    amrex::Array<amrex::Real, NSTATE> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;
    amrex::Array<amrex::Real,NPHASE> AlphaL;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_eL;
    Primitive2State(AlphaL,
                    Alpha_rhoL,
                    Alpha_rho_eL,
       AMREX_D_DECL(uL, vL, wL),
                    pL,
                    PrimitiveL);
    
    amrex::Array<amrex::Real, NSTATE> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_eR;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;
    amrex::Array<amrex::Real,NPHASE> AlphaR;
    Primitive2State(AlphaR,
                    Alpha_rhoR,
                    Alpha_rho_eR,
       AMREX_D_DECL(uR, vR, wR),
                    pR,
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
    
    amrex::Real J_rho_e_mx = 0.0;
    amrex::Array<amrex::Real, NPHASE> J_rho_e;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real Ji = J[i_Phase];
        amrex::Real p_DJi;
        f_Downwind(p_DJi,
                   pL,
                   pR,
                   Ji);
        J_rho_e[i_Phase] = EOS_Parameter.A[i_Phase]*Ji*p_DJi + EOS_Parameter.B[i_Phase]*Ji + EOS_Parameter.D[i_Phase]*J_rho[i_Phase];
        J_rho_e_mx += EOS_Parameter.A[i_Phase]*Ji*p_DJi + EOS_Parameter.B[i_Phase]*Ji + EOS_Parameter.D[i_Phase]*J_rho[i_Phase];
    }
    
    //Flux_PhaseField
    Component2Vector(Flux,
                     J,
                     J_rho,
                     J_rho_e,
        AMREX_D_DECL(M_u, M_v, M_w),
                     M_eK + J_rho_e_mx);
    
}

#endif
#endif