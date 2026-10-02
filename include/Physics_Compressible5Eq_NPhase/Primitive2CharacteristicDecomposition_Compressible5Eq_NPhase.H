// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2CharacteristicDecomposition_Compressible5Eq_NPhase_H_
#define Primitive2CharacteristicDecomposition_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2CharacteristicDecomposition(amrex::Array<amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>,2*NPHASE+AMREX_SPACEDIM+1>& Lambda,
                                           amrex::Array<amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>,2*NPHASE+AMREX_SPACEDIM+1>& Left,
                                           amrex::Array<amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>,2*NPHASE+AMREX_SPACEDIM+1>& Right,
                                           amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Primitive,
                                           Parm const& EOS_Parameter,
                                           int ID_Dimension,
                                           int ID_K)
{
    for (int i_State=0; i_State<2*NPHASE+AMREX_SPACEDIM+1; ++i_State){
        for (int j_State=0; j_State<2*NPHASE+AMREX_SPACEDIM+1; ++j_State){
            Lambda[i_State][j_State] = 0.0;
            Left[i_State][j_State] = 0.0;
            Right[i_State][j_State] = 0.0;
        }
    }
    
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real u;,
                  amrex::Real v;,
	                amrex::Real w );
	  amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Primitive2State(Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha,
                    Primitive);
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    
    amrex::Real c;
    SoundSpeed(c,
               p,
               Alpha,
               Alpha_rho,
               EOS_Parameter,
               ID_K);
    
    amrex::Array<amrex::Real,NPHASE> K;
    NonConservative(K,
                    p,
                    Alpha,
                    Alpha_rho,
                    EOS_Parameter,
                    ID_K);
    amrex::Real Alpha1_rho1 = Alpha_rho[0];
    amrex::Real Alpha2_rho2 = Alpha_rho[1];
    amrex::Real Alpha1 = Alpha[0];
    amrex::Real K1 = K[0];
    // alpha_1 times the Kapila K_1: for NASG with K on, phi_1 (1 + K1) - alpha_1 (see NonConservative)
    amrex::Real Alpha1_K1 = Alpha1*K1;
    if (ID_K == 1 && EOS_Parameter.C[0] > 0.0){
        Alpha1_K1 = NASG_Phi(Alpha1, Alpha1_rho1, EOS_Parameter.C[0])*(1.0 + K1) - Alpha1;
    }
    
    amrex::Real V;
    //x
    if (ID_Dimension == 0){
        V = u;
    }
#if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        V = v;
    }
#endif
#if (AMREX_SPACEDIM > 2)      
    //z
    else if (ID_Dimension == 2){
        V = w;
    }
#endif   
    
    //Lambda
    Lambda[0][0] = V - c;
    Lambda[1][1] = V;
    Lambda[2][2] = V;
#if (AMREX_SPACEDIM > 1)
    Lambda[3][3] = V;
#endif
#if (AMREX_SPACEDIM > 2)    
    Lambda[4][4] = V;
#endif     
    Lambda[2+AMREX_SPACEDIM][2+AMREX_SPACEDIM] = V;
    Lambda[3+AMREX_SPACEDIM][3+AMREX_SPACEDIM] = V + c;
    
    //Left
    Left[0][2+ID_Dimension] = 1.0;
    Left[0][2+AMREX_SPACEDIM] = -1.0/(rho*c);
    
    Left[1][0] = 1.0;
    Left[1][2+AMREX_SPACEDIM] = -Alpha1_rho1/(rho*c*c);
    
    Left[2][1] = 1.0;
    Left[2][2+AMREX_SPACEDIM] = -Alpha2_rho2/(rho*c*c);
    
#if (AMREX_SPACEDIM > 1)
    //x
    if (ID_Dimension == 0){
        Left[3][3] = 1.0;
    }
    //y, z
    else{
        Left[3][2] = 1.0;
    }
#endif
#if (AMREX_SPACEDIM > 2)
    //z
    if (ID_Dimension == 2){
        Left[4][3] = 1.0;
    }
    //x, y
    else{
        Left[4][4] = 1.0;
    }
#endif    
    Left[2+AMREX_SPACEDIM][2+AMREX_SPACEDIM] = (Alpha1_K1)/(rho*c*c);
    Left[2+AMREX_SPACEDIM][3+AMREX_SPACEDIM] = 1.0;
    
    Left[3+AMREX_SPACEDIM][2+ID_Dimension] = 1.0;
    Left[3+AMREX_SPACEDIM][2+AMREX_SPACEDIM] = 1.0/(rho*c);
    
    //Right
    Right[0][0] = -Alpha1_rho1/(2.0*c);
    Right[0][1] = 1.0;
    Right[0][3+AMREX_SPACEDIM] = Alpha1_rho1/(2.0*c);
    
    Right[1][0] = -Alpha2_rho2/(2.0*c);
    Right[1][2] = 1.0;
    Right[1][3+AMREX_SPACEDIM] = Alpha2_rho2/(2.0*c);
    
    //x
    if (ID_Dimension == 0){
        Right[2][0] = 1.0/2.0;
        Right[2][3+AMREX_SPACEDIM] = 1.0/2.0;
    }
    //y, z
    else{
        Right[2][3] = 1.0;
    }
    
#if (AMREX_SPACEDIM > 1)    
    //x
    if (ID_Dimension == 0){
        Right[3][3] = 1.0;
    }
    //y
    else if (ID_Dimension == 1){
        Right[3][0] = 1.0/2.0;
        Right[3][3+AMREX_SPACEDIM] = 1.0/2.0;
    }
    //z
    else if (ID_Dimension == 2){
        Right[3][4] = 1.0;
    }
#endif
#if (AMREX_SPACEDIM > 2)
    //z
    if (ID_Dimension == 2){
        Right[4][0] = 1.0/2.0;
        Right[4][3+AMREX_SPACEDIM] = 1.0/2.0;
    }
    //x, y
    else{
        Right[4][4] = 1.0;
    }
#endif
    Right[2+AMREX_SPACEDIM][0] = -(rho*c)/2.0;
    Right[2+AMREX_SPACEDIM][3+AMREX_SPACEDIM] = (rho*c)/2.0;
    
    Right[3+AMREX_SPACEDIM][0] = (Alpha1_K1)/(2.0*c);
    Right[3+AMREX_SPACEDIM][2+AMREX_SPACEDIM] = 1.0;
    Right[3+AMREX_SPACEDIM][3+AMREX_SPACEDIM] = -(Alpha1_K1)/(2.0*c);
    
}

#endif
#endif