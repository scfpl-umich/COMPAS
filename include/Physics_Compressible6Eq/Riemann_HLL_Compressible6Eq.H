// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLL_Compressible6Eq_H_
#define Riemann_HLL_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

template <typename TF, typename TV, typename TP, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Riemann_HLL(TF& Flux_Riemann,
                 TV& Vel_Riemann,
                 TP& Vel_Alpha1_p1_Riemann,
                 TP& Vel_Alpha2_p2_Riemann,
                 amrex::Real& cmax,
                 TU const& ConservativeL,
                 TU const& ConservativeR,
                 Parm const& EOS_Parameter,
                 int ID_Dimension)
{
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> FluxL;
    Conservative2Flux(FluxL,
                      ConservativeL,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> FluxR;
    Conservative2Flux(FluxR,
                      ConservativeR,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Real Alpha1_rho1L;
    amrex::Real Alpha2_rho2L;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real Alpha1_p1L;
    amrex::Real Alpha2_p2L;
    amrex::Real Alpha1L;
    Primitive2State(Alpha1_rho1L,
                    Alpha2_rho2L,
       AMREX_D_DECL(uL, vL, wL),
                    Alpha1_p1L,
                    Alpha2_p2L,
                    Alpha1L,
                    PrimitiveL);

    amrex::Real Alpha1_rho1R;
    amrex::Real Alpha2_rho2R;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real Alpha1_p1R;
    amrex::Real Alpha2_p2R;
    amrex::Real Alpha1R;
    Primitive2State(Alpha1_rho1R,
                    Alpha2_rho2R,
       AMREX_D_DECL(uR, vR, wR),
                    Alpha1_p1R,
                    Alpha2_p2R,
                    Alpha1R,
                    PrimitiveR);
    
    amrex::Real cL;
    SoundSpeed(cL,
               Alpha1_p1L,
               Alpha2_p2L,
               Alpha1L,
               Alpha1_rho1L,
               Alpha2_rho2L,
               EOS_Parameter);
    
    amrex::Real cR;
    SoundSpeed(cR,
               Alpha1_p1R,
               Alpha2_p2R,
               Alpha1R,
               Alpha1_rho1R,
               Alpha2_rho2R,
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
    
    amrex::Real Sp = std::max( {0.0, VL + cL, VR + cR} );
    amrex::Real Sm = std::min( {0.0, VL - cL, VR - cR} );
    
    Vel_Riemann[0] = (Sp*uL - Sm*uR)/(Sp - Sm);
    Vel_Alpha1_p1_Riemann[0] = (Sp*uL*Alpha1_p1L - Sm*uR*Alpha1_p1R)/(Sp - Sm);
    Vel_Alpha2_p2_Riemann[0] = (Sp*uL*Alpha2_p2L - Sm*uR*Alpha2_p2R)/(Sp - Sm);
#if (AMREX_SPACEDIM > 1)
    Vel_Riemann[1] = (Sp*vL - Sm*vR)/(Sp - Sm);
    Vel_Alpha1_p1_Riemann[1] = (Sp*vL*Alpha1_p1L - Sm*vR*Alpha1_p1R)/(Sp - Sm);
    Vel_Alpha2_p2_Riemann[1] = (Sp*vL*Alpha2_p2L - Sm*vR*Alpha2_p2R)/(Sp - Sm);
#endif
#if (AMREX_SPACEDIM > 2)  
    Vel_Riemann[2] = (Sp*wL - Sm*wR)/(Sp - Sm);
    Vel_Alpha1_p1_Riemann[2] = (Sp*wL*Alpha1_p1L - Sm*wR*Alpha1_p1R)/(Sp - Sm);
    Vel_Alpha2_p2_Riemann[2] = (Sp*wL*Alpha2_p2L - Sm*wR*Alpha2_p2R)/(Sp - Sm);
#endif
    for (int i_State=0; i_State<5+AMREX_SPACEDIM; ++i_State){
        Flux_Riemann[i_State] = ( Sp*FluxL[i_State] - Sm*FluxR[i_State] + Sp*Sm*(ConservativeR[i_State] - ConservativeL[i_State]) )/(Sp - Sm);
    }

    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));

}

#endif
#endif