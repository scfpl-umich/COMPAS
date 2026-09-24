// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLLC_Compressible5Eq_NPhase_H_
#define Riemann_HLLC_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Riemann_HLLC(TF& Flux_Riemann,
                  TV& Vel_Riemann,
                  amrex::Real& cmax,
                  TU const& ConservativeL,
                  TU const& ConservativeR,
                  Parm const& EOS_Parameter,
                  int ID_Dimension,
                  int ID_K)
{

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> FluxL;
    Conservative2Flux(FluxL,
                      ConservativeL,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> FluxR;
    Conservative2Flux(FluxR,
                      ConservativeR,
                      EOS_Parameter,
                      ID_Dimension);
    
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

    // if (Alpha1R >= 0.5){ Alpha1R = 1.0;}
    // else {Alpha1R = 0.0;}
    // if (Alpha1L >= 0.5){ Alpha1L = 1.0;}
    // else {Alpha1L = 0.0;}
    
    amrex::Real cL;
    SoundSpeed(cL,
               pL,
               AlphaL,
               Alpha_rhoL,
               EOS_Parameter,
               ID_K);
    
    amrex::Real cR;
    SoundSpeed(cR,
               pR,
               AlphaR,
               Alpha_rhoR,
               EOS_Parameter,
               ID_K);
    
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
    
    amrex::Real Sm = std::min( 0.0,SL );
    amrex::Real Sp = std::max( 0.0,SR );
    
    amrex::Real rhoL = 0.0;
    amrex::Real rhoR = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rhoL += Alpha_rhoL[i_Phase];
        rhoR += Alpha_rhoR[i_Phase];
    }
    amrex::Real S_Star = ( pR - pL + rhoL*VL*(SL - VL) - rhoR*VR*(SR - VR) )/( rhoL*(SL - VL) - rhoR*(SR - VR) );
    
    amrex::Real sign_S_Star = Sign(S_Star);
    Vel_Riemann[0] = 0.5*(1.0 + sign_S_Star) * uL
                    +0.5*(1.0 - sign_S_Star) * uR;
#if (AMREX_SPACEDIM > 1)    
    Vel_Riemann[1] = 0.5*(1.0 + sign_S_Star) * vL
                    +0.5*(1.0 - sign_S_Star) * vR; 
#endif
#if (AMREX_SPACEDIM > 2)          
    Vel_Riemann[2] = 0.5*(1.0 + sign_S_Star) * wL
                    +0.5*(1.0 - sign_S_Star) * wR;
#endif
    Vel_Riemann[ID_Dimension] += 0.5*(1.0 + sign_S_Star) * Sm * ( (SL - VL)/(SL - S_Star) - 1.0 )
                                +0.5*(1.0 - sign_S_Star) * Sp * ( (SR - VR)/(SR - S_Star) - 1.0 );
    
    amrex::Real ConservativeL_Star;
    amrex::Real ConservativeR_Star;
    for (int i_State=0; i_State<2*NPHASE+AMREX_SPACEDIM+1; ++i_State){
        ConservativeL_Star = (SL - VL)/(SL - S_Star) * ConservativeL[i_State];
        ConservativeR_Star = (SR - VR)/(SR - S_Star) * ConservativeR[i_State];
        if (i_State == NPHASE+ID_Dimension){
            ConservativeL_Star += (SL - VL)/(SL - S_Star) * rhoL * (S_Star - VL);
            ConservativeR_Star += (SR - VR)/(SR - S_Star) * rhoR * (S_Star - VR);
        }
     	  else if (i_State == NPHASE+AMREX_SPACEDIM){
            ConservativeL_Star += (SL - VL)/(SL - S_Star) * (S_Star - VL) * ( rhoL*S_Star + pL/(SL - VL) );
            ConservativeR_Star += (SR - VR)/(SR - S_Star) * (S_Star - VR) * ( rhoR*S_Star + pR/(SR - VR) );
        }
        Flux_Riemann[i_State] = 0.5*(1.0 + sign_S_Star) * ( FluxL[i_State] + Sm*(ConservativeL_Star - ConservativeL[i_State]) )
                               +0.5*(1.0 - sign_S_Star) * ( FluxR[i_State] + Sp*(ConservativeR_Star - ConservativeR[i_State]) );
    }
    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));
}

#endif
#endif