// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_LLF_Compressible5Eq_NPhase_H_
#define Riemann_LLF_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Riemann_LLF(TF& Flux_Riemann,
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
    
    amrex::Real Lambda = std::max( std::abs(VL) + cL, std::abs(VR) + cR );
    
    Vel_Riemann[0] = (uL + uR)/2.0;
#if (AMREX_SPACEDIM > 1)
    Vel_Riemann[1] = (vL + vR)/2.0;
#endif
#if (AMREX_SPACEDIM > 2)  
    Vel_Riemann[2] = (wL + wR)/2.0;
#endif
    for (int i_State=0; i_State<2*NPHASE+AMREX_SPACEDIM+1; ++i_State){
        Flux_Riemann[i_State] = ( FluxL[i_State] + FluxR[i_State] - Lambda*(ConservativeR[i_State] - ConservativeL[i_State]) )/2.0;
    }

    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));

}

#endif
#endif