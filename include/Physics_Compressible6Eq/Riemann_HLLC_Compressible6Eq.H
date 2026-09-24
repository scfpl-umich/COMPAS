// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Riemann_HLLC_Compressible6Eq_H_
#define Riemann_HLLC_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

template <typename TF, typename TV, typename TP, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Riemann_HLLC(TF& Flux_Riemann,
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
    AMREX_D_TERM(amrex::Real uL;,
                 amrex::Real vL;,
                 amrex::Real wL);
    amrex::Real Alpha1_p1L;
    amrex::Real Alpha2_p2L;
    amrex::Real Alpha1L;
    amrex::Real Alpha2L;
    Primitive2State(Alpha1_rho1L,
                    Alpha2_rho2L,
       AMREX_D_DECL(uL, vL, wL),
                    Alpha1_p1L,
                    Alpha2_p2L,
                    Alpha1L,
                    PrimitiveL);
    Alpha2L = 1.0 - Alpha1L;

    amrex::Real Alpha1_rho1R;
    amrex::Real Alpha2_rho2R;
    AMREX_D_TERM(amrex::Real uR;,
                 amrex::Real vR;,
                 amrex::Real wR);
    amrex::Real Alpha1_p1R;
    amrex::Real Alpha2_p2R;
    amrex::Real Alpha1R;
    amrex::Real Alpha2R;
    Primitive2State(Alpha1_rho1R,
                    Alpha2_rho2R,
       AMREX_D_DECL(uR, vR, wR),
                    Alpha1_p1R,
                    Alpha2_p2R,
                    Alpha1R,
                    PrimitiveR);
    Alpha2R = 1.0 - Alpha1R;

    AMREX_D_TERM(amrex::Real rho_uL;,
                 amrex::Real rho_vL;,
                 amrex::Real rho_wL);
    amrex::Real Alpha1_rho1_E1L;
    amrex::Real Alpha2_rho2_E2L;
    Conservative2State(Alpha1_rho1L,
                       Alpha2_rho2L,
          AMREX_D_DECL(rho_uL, rho_vL, rho_wL),
                       Alpha1_rho1_E1L,
                       Alpha2_rho2_E2L,
                       Alpha1L,
                       ConservativeL);

    AMREX_D_TERM(amrex::Real rho_uR;,
                 amrex::Real rho_vR;,
                 amrex::Real rho_wR);
    amrex::Real Alpha1_rho1_E1R;
    amrex::Real Alpha2_rho2_E2R;
    Conservative2State(Alpha1_rho1R,
                       Alpha2_rho2R,
          AMREX_D_DECL(rho_uR, rho_vR, rho_wR),
                       Alpha1_rho1_E1R,
                       Alpha2_rho2_E2R,
                       Alpha1R,
                       ConservativeR);

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

    amrex::Real SL = std::min(VL - cL, VR - cR);
    amrex::Real SR = std::max(VL + cL, VR + cR);
    amrex::Real Sm = std::min(0.0, SL);
    amrex::Real Sp = std::max(0.0, SR);

    amrex::Real rhoL = Alpha1_rho1L
                     + Alpha2_rho2L;
    amrex::Real rhoR = Alpha1_rho1R
                     + Alpha2_rho2R;
    amrex::Real pL   = Alpha1_p1L
                     + Alpha2_p2L;
    amrex::Real pR   = Alpha1_p1R
                     + Alpha2_p2R;
    amrex::Real S_Star = (pR - pL + rhoL*VL*(SL - VL) - rhoR*VR*(SR - VR))/(rhoL*(SL - VL) - rhoR*(SR - VR));

    amrex::Real Alpha1_rho1L_Star = (SL - VL)/(SL - S_Star)*Alpha1_rho1L;
    amrex::Real Alpha1_rho1R_Star = (SR - VR)/(SR - S_Star)*Alpha1_rho1R;
    amrex::Real Alpha2_rho2L_Star = (SL - VL)/(SL - S_Star)*Alpha2_rho2L;
    amrex::Real Alpha2_rho2R_Star = (SR - VR)/(SR - S_Star)*Alpha2_rho2R;

    amrex::Real rho_uL_Star = (SL - VL)/(SL - S_Star)*rhoL*(uL + amrex::Real(ID_Dimension == 0)*(S_Star - uL));
    amrex::Real rho_uR_Star = (SR - VR)/(SR - S_Star)*rhoR*(uR + amrex::Real(ID_Dimension == 0)*(S_Star - uR));
#if (AMREX_SPACEDIM > 1)
    amrex::Real rho_vL_Star = (SL - VL)/(SL - S_Star)*rhoL*(vL + amrex::Real(ID_Dimension == 1)*(S_Star - vL));
    amrex::Real rho_vR_Star = (SR - VR)/(SR - S_Star)*rhoR*(vR + amrex::Real(ID_Dimension == 1)*(S_Star - vR));
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real rho_wL_Star = (SL - VL)/(SL - S_Star)*rhoL*(wL + amrex::Real(ID_Dimension == 2)*(S_Star - wL));
    amrex::Real rho_wR_Star = (SR - VR)/(SR - S_Star)*rhoR*(wR + amrex::Real(ID_Dimension == 2)*(S_Star - wR));
#endif

    amrex::Real Alpha1_p1L_Star = Alpha1_p1L + (VL - SL)*Alpha1_rho1L*VL - (S_Star - SL)*Alpha1_rho1L_Star*S_Star;
    amrex::Real Alpha1_p1R_Star = Alpha1_p1R + (VR - SR)*Alpha1_rho1R*VR - (S_Star - SR)*Alpha1_rho1R_Star*S_Star;
    amrex::Real Alpha2_p2L_Star = Alpha2_p2L + (VL - SL)*Alpha2_rho2L*VL - (S_Star - SL)*Alpha2_rho2L_Star*S_Star;
    amrex::Real Alpha2_p2R_Star = Alpha2_p2R + (VR - SR)*Alpha2_rho2R*VR - (S_Star - SR)*Alpha2_rho2R_Star*S_Star;

    amrex::Real Alpha1_rho1_E1L_Star = ((VL - SL)*Alpha1_rho1_E1L + VL*Alpha1_p1L - S_Star*Alpha1_p1L_Star)/(S_Star - SL);
    amrex::Real Alpha1_rho1_E1R_Star = ((VR - SR)*Alpha1_rho1_E1R + VR*Alpha1_p1R - S_Star*Alpha1_p1R_Star)/(S_Star - SR);
    amrex::Real Alpha2_rho2_E2L_Star = ((VL - SL)*Alpha2_rho2_E2L + VL*Alpha2_p2L - S_Star*Alpha2_p2L_Star)/(S_Star - SL);
    amrex::Real Alpha2_rho2_E2R_Star = ((VR - SR)*Alpha2_rho2_E2R + VR*Alpha2_p2R - S_Star*Alpha2_p2R_Star)/(S_Star - SR);

    amrex::Real Alpha1L_Star = Alpha1L;
    amrex::Real Alpha1R_Star = Alpha1R;

    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> ConservativeL_Star;
    State2Conservative(ConservativeL_Star,
                       Alpha1_rho1L_Star,
                       Alpha2_rho2L_Star,
          AMREX_D_DECL(rho_uL_Star, rho_vL_Star, rho_wL_Star),
                       Alpha1_rho1_E1L_Star,
                       Alpha2_rho2_E2L_Star,
                       Alpha1L_Star);

    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> ConservativeR_Star;
    State2Conservative(ConservativeR_Star,
                       Alpha1_rho1R_Star,
                       Alpha2_rho2R_Star,
          AMREX_D_DECL(rho_uR_Star, rho_vR_Star, rho_wR_Star),
                       Alpha1_rho1_E1R_Star,
                       Alpha2_rho2_E2R_Star,
                       Alpha1R_Star);

    amrex::Real sign_S_Star = Sign(S_Star);
    Vel_Riemann[0] = 0.5*(1.0 + sign_S_Star)*uL
                   + 0.5*(1.0 - sign_S_Star)*uR;
    Vel_Alpha1_p1_Riemann[0] = 0.5*(1.0 + sign_S_Star)*(uL*Alpha1_p1L + amrex::Real(ID_Dimension == 0)*Sm*(VL*Alpha1_p1L - S_Star*Alpha1_p1L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(uR*Alpha1_p1R + amrex::Real(ID_Dimension == 0)*Sp*(VR*Alpha1_p1R - S_Star*Alpha1_p1R_Star)/(S_Star - SR));
    Vel_Alpha2_p2_Riemann[0] = 0.5*(1.0 + sign_S_Star)*(uL*Alpha2_p2L + amrex::Real(ID_Dimension == 0)*Sm*(VL*Alpha2_p2L - S_Star*Alpha2_p2L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(uR*Alpha2_p2R + amrex::Real(ID_Dimension == 0)*Sp*(VR*Alpha2_p2R - S_Star*Alpha2_p2R_Star)/(S_Star - SR));
#if (AMREX_SPACEDIM > 1)
    Vel_Riemann[1] = 0.5*(1.0 + sign_S_Star)*vL
                   + 0.5*(1.0 - sign_S_Star)*vR;
    Vel_Alpha1_p1_Riemann[1] = 0.5*(1.0 + sign_S_Star)*(vL*Alpha1_p1L + amrex::Real(ID_Dimension == 1)*Sm*(VL*Alpha1_p1L - S_Star*Alpha1_p1L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(vR*Alpha1_p1R + amrex::Real(ID_Dimension == 1)*Sp*(VR*Alpha1_p1R - S_Star*Alpha1_p1R_Star)/(S_Star - SR));
    Vel_Alpha2_p2_Riemann[1] = 0.5*(1.0 + sign_S_Star)*(vL*Alpha2_p2L + amrex::Real(ID_Dimension == 1)*Sm*(VL*Alpha2_p2L - S_Star*Alpha2_p2L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(vR*Alpha2_p2R + amrex::Real(ID_Dimension == 1)*Sp*(VR*Alpha2_p2R - S_Star*Alpha2_p2R_Star)/(S_Star - SR));
#endif
#if (AMREX_SPACEDIM > 2)
    Vel_Riemann[2] = 0.5*(1.0 + sign_S_Star)*wL
                   + 0.5*(1.0 - sign_S_Star)*wR;
    Vel_Alpha1_p1_Riemann[2] = 0.5*(1.0 + sign_S_Star)*(wL*Alpha1_p1L + amrex::Real(ID_Dimension == 2)*Sm*(VL*Alpha1_p1L - S_Star*Alpha1_p1L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(wR*Alpha1_p1R + amrex::Real(ID_Dimension == 2)*Sp*(VR*Alpha1_p1R - S_Star*Alpha1_p1R_Star)/(S_Star - SR));
    Vel_Alpha2_p2_Riemann[2] = 0.5*(1.0 + sign_S_Star)*(wL*Alpha2_p2L + amrex::Real(ID_Dimension == 2)*Sm*(VL*Alpha2_p2L - S_Star*Alpha2_p2L_Star)/(S_Star - SL))
                             + 0.5*(1.0 - sign_S_Star)*(wR*Alpha2_p2R + amrex::Real(ID_Dimension == 2)*Sp*(VR*Alpha2_p2R - S_Star*Alpha2_p2R_Star)/(S_Star - SR));
#endif

    for (int i_State = 0; i_State < 5+AMREX_SPACEDIM; ++i_State){
        Flux_Riemann[i_State] = 0.5*(1.0 + sign_S_Star)*(FluxL[i_State] + Sm*(ConservativeL_Star[i_State] - ConservativeL[i_State]))
                              + 0.5*(1.0 - sign_S_Star)*(FluxR[i_State] + Sp*(ConservativeR_Star[i_State] - ConservativeR[i_State]));
    }

    cmax = std::max(cL + std::sqrt(AMREX_D_TERM(uL*uL, + vL*vL, + wL*wL)), cR + std::sqrt(AMREX_D_TERM(uR*uR, + vR*vR, + wR*wR)));

}

#endif
#endif