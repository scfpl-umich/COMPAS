// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_LLF_Compressible6Eq_IE_NPhase_H_
#define Riemann_LLF_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
Riemann_LLF (TF & Flux_Riemann,
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
    Conservative2State(AlphaL, Alpha_rhoL, Alpha_rho_eL, AMREX_D_DECL(rho_uL, rho_vL, rho_wL), rho_EL, ConservativeL);

    
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

    Vel_Riemann[0] = 0.5* (uL+uR);
    #if (AMREX_SPACEDIM > 1)    
        Vel_Riemann[1] = 0.5* (vL+vR);
    #endif
    #if (AMREX_SPACEDIM > 2)          
        Vel_Riemann[2] = 0.5* (wL+wR);
    #endif

    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));


    for (int i_State=0; i_State < NSTATE; ++i_State){

        Flux_Riemann[i_State] = 0.5*( FluxL[i_State] + FluxR[i_State]) - 0.5*cmax*(ConservativeR[i_State] - ConservativeL[i_State]);
    }
    
}

#endif
#endif