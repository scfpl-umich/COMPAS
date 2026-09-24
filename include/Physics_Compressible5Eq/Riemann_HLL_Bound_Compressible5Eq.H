// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLL_Bound_Compressible5Eq_H_
#define Riemann_HLL_Bound_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Riemann_HLL_Bound(TF& Flux_Riemann,
                       TV& Vel_Riemann,
                       amrex::Real& cmax,
                       TU const& ConservativeL,
                       TU const& ConservativeR,
                       TU const& ConservativeL_c,
                       TU const& ConservativeR_c,
                       Parm const& EOS_Parameter,
                       int ID_Dimension,
                       int ID_K,
                       int ID_Lambda)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> FluxL;
    Conservative2Flux(FluxL,
                      ConservativeL,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> FluxR;
    Conservative2Flux(FluxR,
                      ConservativeR,
                      EOS_Parameter,
                      ID_Dimension);
    
    amrex::Real Alpha1_rho1L;
    amrex::Real Alpha2_rho2L;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;
    amrex::Real Alpha1L;
    Primitive2State(Alpha1_rho1L,
                    Alpha2_rho2L,
       AMREX_D_DECL(uL, vL, wL),
                    pL,
                    Alpha1L,
                    PrimitiveL);

    amrex::Real Alpha1_rho1R;
    amrex::Real Alpha2_rho2R;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;
    amrex::Real Alpha1R;
    Primitive2State(Alpha1_rho1R,
                    Alpha2_rho2R,
       AMREX_D_DECL(uR, vR, wR),
                    pR,
                    Alpha1R,
                    PrimitiveR);
    
    amrex::Real LambdaLL_p;
    amrex::Real LambdaLL_m;
    Conservative2Lambda_HLL(LambdaLL_p, LambdaLL_m,
                            ConservativeL,
                            ConservativeL_c,
                            EOS_Parameter,
                            ID_Dimension,
                            ID_K,
                            ID_Lambda);
    amrex::Real LambdaLR_p;
    amrex::Real LambdaLR_m;
    Conservative2Lambda_HLL(LambdaLR_p, LambdaLR_m,
                            ConservativeL,
                            ConservativeR_c,
                            EOS_Parameter,
                            ID_Dimension,
                            ID_K,
                            ID_Lambda);
    amrex::Real LambdaL_p = std::max(LambdaLL_p, LambdaLR_p);
    amrex::Real LambdaL_m = std::min(LambdaLL_m, LambdaLR_m);
    
    amrex::Real LambdaRL_p;
    amrex::Real LambdaRL_m;
    Conservative2Lambda_HLL(LambdaRL_p, LambdaRL_m,
                            ConservativeR,
                            ConservativeL_c,
                            EOS_Parameter,
                            ID_Dimension,
                            ID_K,
                            ID_Lambda);
    amrex::Real LambdaRR_p;
    amrex::Real LambdaRR_m;
    Conservative2Lambda_HLL(LambdaRR_p, LambdaRR_m,
                            ConservativeR,
                            ConservativeR_c,
                            EOS_Parameter,
                            ID_Dimension,
                            ID_K,
                            ID_Lambda);
    amrex::Real LambdaR_p = std::max(LambdaRL_p, LambdaRR_p);
    amrex::Real LambdaR_m = std::min(LambdaRL_m, LambdaRR_m);
    
    amrex::Real Sp = std::max({ 0.0, LambdaL_p, LambdaR_p }) + 1.0e-10;
    amrex::Real Sm = std::min({ 0.0, LambdaL_m, LambdaR_m }) - 1.0e-10;
    
    Vel_Riemann[0] = (Sp*uL - Sm*uR)/(Sp - Sm);
#if (AMREX_SPACEDIM > 1)
    Vel_Riemann[1] = (Sp*vL - Sm*vR)/(Sp - Sm);
#endif
#if (AMREX_SPACEDIM > 2)  
    Vel_Riemann[2] = (Sp*wL - Sm*wR)/(Sp - Sm);
#endif
    for (int i_State=0; i_State<4+AMREX_SPACEDIM; ++i_State){
        Flux_Riemann[i_State] = ( Sp*FluxL[i_State] - Sm*FluxR[i_State] + Sp*Sm*(ConservativeR[i_State] - ConservativeL[i_State]) )/(Sp - Sm);
    }

    cmax = std::max( -Sm, Sp );

}

#endif
#endif