// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLL_Compressible6Eq_IE_NPhase_H_
#define Riemann_HLL_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
Riemann_HLL (TF & Flux_Riemann,
             TV & Vel_Riemann,
             amrex::Real & cmax,
             TU const & ConservativeL,
             TU const & ConservativeR,
             Parm const& EOS_Parameter,
             int ID_Dimension)
{
    amrex::Array<amrex::Real,NSTATE> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,NSTATE> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,NSTATE> FluxL;
    Conservative2Flux(FluxL,
                      ConservativeL,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real,NSTATE> FluxR;
    Conservative2Flux(FluxR,
                      ConservativeR,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real,NPHASE> AlphaL;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL;
    amrex::Array<amrex::Real,NPHASE> Alpha_PsL;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;

    Primitive2State(AlphaL, Alpha_rhoL, Alpha_PsL,
       AMREX_D_DECL(uL, vL, wL),
                    pL,
                    PrimitiveL);

    amrex::Array<amrex::Real,NPHASE> AlphaR;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR;
    amrex::Array<amrex::Real,NPHASE> Alpha_PsR;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;
 
    Primitive2State(AlphaR, Alpha_rhoR, Alpha_PsR,
       AMREX_D_DECL(uR, vR, wR),
                    pR,
                    PrimitiveR);
    
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
    
    amrex::Real Sp = std::max(std::max(0.0, VL + cL), VR + cR );
    amrex::Real Sm = std::min(std::min(0.0, VL - cL), VR - cR );

    // amrex::Real Sp = std::max( {VL,VR,VL - cL, VR - cR, VL + cL, VR + cR} );


    Vel_Riemann[0] = (Sp*uL - Sm*uR)/(Sp - Sm);
#if (AMREX_SPACEDIM > 1)
    Vel_Riemann[1] = (Sp*vL - Sm*vR)/(Sp - Sm);
#endif
#if (AMREX_SPACEDIM > 2)  
    Vel_Riemann[2] = (Sp*wL - Sm*wR)/(Sp - Sm);
#endif
    for (int i_State=0; i_State< NSTATE; ++i_State){
        Flux_Riemann[i_State] = ( Sp*FluxL[i_State] - Sm*FluxR[i_State] + Sp*Sm*(ConservativeR[i_State] - ConservativeL[i_State]) )/(Sp - Sm);
    }

    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));
}

#endif
#endif