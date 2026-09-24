// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLLC_Compressible6Eq_IE_NPhase_H_
#define Riemann_HLLC_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
Riemann_HLLC (TF & Flux_Riemann,
              TV & Vel_Riemann,
              amrex::Real & cmax,
              TU const & ConservativeL,
              TU const & ConservativeR,
              Parm const& EOS_Parameter,
              int ID_Dimension)
{

    amrex::Array<amrex::Real, NSTATE> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real, NSTATE> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real, NSTATE> FluxL;
    Conservative2Flux(FluxL,
                      ConservativeL,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real, NSTATE> FluxR;
    Conservative2Flux(FluxR,
                      ConservativeR,
                      EOS_Parameter,
                      ID_Dimension);


    // data on the left
    amrex::Array<amrex::Real,NPHASE> AlphaL;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL;
    amrex::Array<amrex::Real,NPHASE> Alpha_PsL;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_eL;

    AMREX_D_TERM( amrex::Real rho_uL;, 
                  amrex::Real rho_vL;, 
                  amrex::Real rho_wL);
                
    amrex::Real rho_EL;

    Primitive2State(AlphaL, Alpha_rhoL, Alpha_PsL, AMREX_D_DECL(uL, vL, wL), pL, PrimitiveL);
    Conservative2State(AlphaL, Alpha_rhoL, Alpha_rho_eL, AMREX_D_DECL(rho_uL, rho_vL, rho_wL),rho_EL, ConservativeL);

    
    // data on the right
    amrex::Array<amrex::Real,NPHASE> AlphaR;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR;
    amrex::Array<amrex::Real,NPHASE> Alpha_PsR;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_eR;

    AMREX_D_TERM( amrex::Real rho_uR;, 
                  amrex::Real rho_vR;, 
                  amrex::Real rho_wR);

    amrex::Real rho_ER;

    Primitive2State(AlphaR, Alpha_rhoR, Alpha_PsR, AMREX_D_DECL(uR, vR, wR), pR, PrimitiveR);
    Conservative2State(AlphaR, Alpha_rhoR, Alpha_rho_eR, AMREX_D_DECL(rho_uR, rho_vR, rho_wR),rho_ER, ConservativeR);

    amrex::Real cL;
    SoundSpeed(cL,
                AlphaL,
                Alpha_rhoL,
                Alpha_PsL,
                EOS_Parameter);
    
    amrex::Real cR;
    SoundSpeed(cR,
        AlphaR,
        Alpha_rhoR,
        Alpha_PsR,
        EOS_Parameter);


    amrex::Real VL;
    amrex::Real VR;
	//x
    if (ID_Dimension == 0){
        VL = uL;
        VR = uR;
    }
    #if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        VL = vL;
        VR = vR;
    }
    #endif
    #if (AMREX_SPACEDIM > 2)      
    //z
    else if (ID_Dimension == 2){
        VL = wL;
        VR = wR;
    }
    #endif
    
    amrex::Real V_avg = 0.5*(VL + VR);
    amrex::Real c_avg = 0.5*(cL + cR);
    amrex::Real SL = std::min( V_avg-c_avg,VL-cL );
    amrex::Real SR = std::max( V_avg+c_avg,VR+cR );
    
    amrex::Real Sm = std::min( 0.0, SL );
    amrex::Real Sp = std::max( 0.0, SR );

    amrex::Array<amrex::Real, NSTATE> QR, QL, ConservativeL_Star, ConservativeR_Star;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> VelL = {AMREX_D_DECL(uL, vL, wL)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> VelR = {AMREX_D_DECL(uR, vR, wR)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> VelL_Star, VelR_Star;
    amrex::Real S_Star, P_Star;
    

    for(int i=0; i<NSTATE; ++i){
        QR[i] = FluxR[i] - SR*ConservativeR[i];
        QL[i] = FluxL[i] - SL*ConservativeL[i];
        ConservativeL_Star[i] = ConservativeL[i];
        ConservativeR_Star[i] = ConservativeR[i];

    }


    amrex::Real QL_rho = 0.;
    amrex::Real QR_rho = 0.;
    
    for(int i_Phase=0; i_Phase < NPHASE; ++i_Phase){
        QL_rho += QL[NPHASE + i_Phase];
        QR_rho += QR[NPHASE + i_Phase];
    }
    

    for(int i=0; i<AMREX_SPACEDIM; ++i){
        VelL_Star[i] = VelL[i];
        VelR_Star[i] = VelR[i];
    }

    S_Star = (QL[3*NPHASE+ ID_Dimension] - QR[3*NPHASE+ID_Dimension])/(QL_rho - QR_rho);
    P_Star = -(QL[3*NPHASE+ ID_Dimension]*QR_rho - QR[3*NPHASE+ID_Dimension]*QL_rho)/(QL_rho - QR_rho);

    VelL_Star[ID_Dimension] = S_Star;
    VelR_Star[ID_Dimension] = S_Star;
    
    amrex::Real sign_S_Star = Sign(S_Star);

    for(int id = 0; id<AMREX_SPACEDIM; ++id){
        Vel_Riemann[id] = 0.5*(1.0 + sign_S_Star) *( VelL[id] + Sm * ( (SL - VelL[id])/(SL - VelL_Star[id]) - 1.0 ))
                    +0.5*(1.0 - sign_S_Star) * ( VelR[id] + Sp * ( (SR - VelR[id])/(SR - VelR_Star[id]) - 1.0 ));
    }

    // Vel_Riemann[ID_Dimension] = 0.5*(1.0 + sign_S_Star) *( VL + Sm * ( (SL - VL)/(SL - S_Star) - 1.0 ))
    //                 +0.5*(1.0 - sign_S_Star) * ( VR + Sp * ( (SR - VR)/(SR - S_Star) - 1.0 ));

    amrex::Real XiL, XiR;

    // amrex::Real S_Star = VelL_Star[ID_Dimension];
    
    // amrex::Real sign_S_Star = Sign(S_Star);

    XiL = (VelL[ID_Dimension] - SL)/(S_Star - SL);
    XiR = (VelR[ID_Dimension] - SR)/(S_Star - SR);

    // alpha_k
    amrex::Array<amrex::Real,NPHASE> AlphaL_Star, AlphaR_Star;
    for(int i_state=0; i_state<NPHASE; ++i_state){
        ConservativeL_Star[i_state] = QL[i_state]/(S_Star - SL);
        ConservativeR_Star[i_state] = QR[i_state]/(S_Star - SR);
        
        AlphaL_Star[i_state] = ConservativeL_Star[i_state];
        AlphaR_Star[i_state] = ConservativeR_Star[i_state];
    }

    // alpha_k rho_k
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL_Star, Alpha_rhoR_Star;
    for(int i_state=NPHASE; i_state<2*NPHASE; ++i_state){
        ConservativeL_Star[i_state] = QL[i_state]/(S_Star - SL);
        ConservativeR_Star[i_state] = QR[i_state]/(S_Star - SR);
        int ik = i_state - NPHASE;
        Alpha_rhoL_Star[ik] = ConservativeL_Star[i_state];
        Alpha_rhoR_Star[ik] = ConservativeR_Star[i_state];
    }

    // alpha_k rho_k e_k
    amrex::Array<amrex::Real,NPHASE> Alpha_PsL_Star, Alpha_PsR_Star;
    for(int i_state=2*NPHASE; i_state<3*NPHASE; ++i_state){
        
        int ik = i_state - 2*NPHASE;
        
        amrex::Real A = EOS_Parameter.A[ik];
        amrex::Real B = EOS_Parameter.B[ik];

        // Alpha_PsL_Star[ik] = (A*Alpha_PsL[ik] + AlphaL[ik]*B - AlphaL_Star[ik]*B*(Alpha_rhoL[ik]/Alpha_rhoL_Star[ik]) 
        //                     - 0.5*Alpha_PsL[ik]*(Alpha_rhoL[ik]/Alpha_rhoL_Star[ik] - 1.))/(A + 0.5*(Alpha_rhoL[ik]/Alpha_rhoL_Star[ik] -1.));

        // Alpha_PsR_Star[ik] = (A*Alpha_PsR[ik] + AlphaR[ik]*B - AlphaR_Star[ik]*B*(Alpha_rhoR[ik]/Alpha_rhoR_Star[ik]) 
        //                     - 0.5*Alpha_PsR[ik]*(Alpha_rhoR[ik]/Alpha_rhoR_Star[ik] - 1.))/(A + 0.5*(Alpha_rhoR[ik]/Alpha_rhoR_Star[ik] -1.));

        Alpha_PsL_Star[ik] = (A*Alpha_PsL[ik] + AlphaL[ik]*B - AlphaL_Star[ik]*B*(1./XiL) 
                            - 0.5*Alpha_PsL[ik]*(1./XiL - 1.))/(A + 0.5*(1./XiL -1.));

        Alpha_PsR_Star[ik] = (A*Alpha_PsR[ik] + AlphaR[ik]*B - AlphaR_Star[ik]*B*(1./XiR) 
                            - 0.5*Alpha_PsR[ik]*(1./XiR - 1.))/(A + 0.5*(1./XiR -1.));
        

        bool PressureBase = true;//Pressure
        bool OutputInternalEnergy = true;//From Pressure to rho_e 
        amrex::Real akrhokekL_star, akrhokekR_star;

        EquationOfState_Phase(akrhokekL_star,
                    Alpha_PsL_Star[ik],
                    AlphaL_Star[ik],
                    Alpha_rhoL_Star[ik],
                    EOS_Parameter,
                    ik,
                    PressureBase,
                    OutputInternalEnergy);

        EquationOfState_Phase(akrhokekR_star,
                Alpha_PsR_Star[ik],
                AlphaR_Star[ik],
                Alpha_rhoR_Star[ik],
                EOS_Parameter,
                ik,
                PressureBase,
                OutputInternalEnergy);
        
    
        ConservativeL_Star[i_state] = akrhokekL_star;
        ConservativeR_Star[i_state] = akrhokekR_star;
    }

    // rho u
    for(int i_state=3*NPHASE; i_state<3*NPHASE + AMREX_SPACEDIM; ++i_state){
        ConservativeL_Star[i_state] = VelL_Star[i_state - 3*NPHASE]*QL_rho/(S_Star - SL);
        ConservativeR_Star[i_state] = VelR_Star[i_state - 3*NPHASE]*QR_rho/(S_Star - SR);
    }

    // rho E
    ConservativeL_Star[3*NPHASE + AMREX_SPACEDIM] = (QL[3*NPHASE + AMREX_SPACEDIM] - P_Star*S_Star)/(S_Star - SL);
    ConservativeR_Star[3*NPHASE + AMREX_SPACEDIM] = (QR[3*NPHASE + AMREX_SPACEDIM] - P_Star*S_Star)/(S_Star - SR);


    for (int i_State=0; i_State < NSTATE; ++i_State){
        Flux_Riemann[i_State] = 0.5*(1.0 + sign_S_Star) * ( FluxL[i_State] + Sm*(ConservativeL_Star[i_State] - ConservativeL[i_State]) )
                               +0.5*(1.0 - sign_S_Star) * ( FluxR[i_State] + Sp*(ConservativeR_Star[i_State] - ConservativeR[i_State]) );

    }
    
    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));
   
}

#endif
#endif