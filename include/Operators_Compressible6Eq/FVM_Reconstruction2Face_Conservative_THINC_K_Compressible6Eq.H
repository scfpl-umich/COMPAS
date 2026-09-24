// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_THINC_K_Compressible6Eq_H_
#define FVM_Reconstruction2Face_Conservative_THINC_K_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
Reconstruct_State_From_VolumeFraction (amrex::Array<amrex::Real,5+AMREX_SPACEDIM> & Conservative_Reconstructed,
                                       amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const & Conservative,
                                       amrex::Real Alpha1_Reconstructed,
                                       Parm const& EOS_Parameter)
{

    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;)
    amrex::Real Alpha1_p1;
    amrex::Real Alpha2_p2;
    amrex::Real Alpha1;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive;

    Conservative2Primitive(Primitive,Conservative,EOS_Parameter);
    Primitive2State(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u,v,w), Alpha1_p1, Alpha2_p2, Alpha1, Primitive);
    
    amrex::Real rho1 = Alpha1_rho1 / Alpha1;
    amrex::Real rho2 = Alpha2_rho2 / (1.0 - Alpha1);
    amrex::Real Alpha1_rho1_Recon = Alpha1_Reconstructed * rho1;
    amrex::Real Alpha2_rho2_Recon = (1.0 - Alpha1_Reconstructed) * rho2;
    amrex::Real rho_Recon = Alpha1_rho1_Recon + Alpha2_rho2_Recon;
    
    AMREX_D_TERM(
    amrex::Real Mom_x_Recon = u * rho_Recon;,
    amrex::Real Mom_y_Recon = v * rho_Recon;,
    amrex::Real Mom_z_Recon = w * rho_Recon;)

    amrex::Real p1 = Alpha1_p1 / Alpha1;
    amrex::Real p2 = Alpha2_p2 / (1.0 - Alpha1);

    amrex::Real K = 0.5 * (AMREX_D_TERM(u*u, + v*v, + w*w));
    amrex::Real rho1_e1 = EOS_Parameter.A1 * p1 + EOS_Parameter.B1 + rho1 * EOS_Parameter.D1;
    amrex::Real rho2_e2 = EOS_Parameter.A2 * p2 + EOS_Parameter.B2 + rho2 * EOS_Parameter.D2;

    Conservative_Reconstructed[INDEX_Mass1] = Alpha1_rho1_Recon;
    Conservative_Reconstructed[INDEX_Mass2] = Alpha2_rho2_Recon;
    AMREX_D_TERM(
    Conservative_Reconstructed[INDEX_MomentumX] = u * rho_Recon;,
    Conservative_Reconstructed[INDEX_MomentumY] = v * rho_Recon;,
    Conservative_Reconstructed[INDEX_MomentumZ] = w * rho_Recon;)
    Conservative_Reconstructed[INDEX_Energy1] = K * Alpha1_rho1_Recon + rho1_e1 * Alpha1_Reconstructed;
    Conservative_Reconstructed[INDEX_Energy2] = K * Alpha2_rho2_Recon + rho2_e2 * (1.0 - Alpha1_Reconstructed);
    Conservative_Reconstructed[INDEX_VolumeFraction1] = Alpha1_Reconstructed;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
FVM_Reconstruction2Face_Conservative_THINC_K(int i, int j, int k,
                                             amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeL,
                                             amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                             amrex::Array4<amrex::Real const> const& Conservative,
                                             Parm const& EOS_Parameter,
                                             amrex::Array<amrex::Real, AMREX_SPACEDIM> const& dx,
                                             int ID_Dimension,
                                             int ID_THINC_UseNormal,
                                             amrex::Real NormalThreshold,
                                             amrex::Real EPS, 
                                             amrex::Real Beta,
                                             int ID_WENO, int ID_Limiter, int ID_Poly)
{


    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_mm;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_m;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_c;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_p;

    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal;

    for (int i_State = 0; i_State < 5+AMREX_SPACEDIM; i_State++){
        Conservative_c[i_State] = Conservative(i,j,k,i_State);
        //x
        if (ID_Dimension == 0){
            Conservative_mm[i_State] = Conservative(i-2,j,k,i_State);
            Conservative_m [i_State] = Conservative(i-1,j,k,i_State);
            Conservative_p [i_State] = Conservative(i+1,j,k,i_State);
            FDM_Normal_XFaceFromCenter_K(i,j,k,INDEX_VolumeFraction1, Normal, Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));

        }
        //y
        else if (ID_Dimension == 1){
            Conservative_mm[i_State] = Conservative(i,j-2,k,i_State);
            Conservative_m [i_State] = Conservative(i,j-1,k,i_State);
            Conservative_p [i_State] = Conservative(i,j+1,k,i_State);
#if (AMREX_SPACEDIM > 1)
            FDM_Normal_YFaceFromCenter_K(i,j,k,INDEX_VolumeFraction1, Normal, Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
        //z
        else if (ID_Dimension == 2){
            Conservative_mm[i_State] = Conservative(i,j,k-2,i_State);
            Conservative_m [i_State] = Conservative(i,j,k-1,i_State);
            Conservative_p [i_State] = Conservative(i,j,k+1,i_State);
#if (AMREX_SPACEDIM > 2)
            FDM_Normal_ZFaceFromCenter_K(i,j,k,INDEX_VolumeFraction1, Normal, Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
    }

    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive_mm, Primitive_m, Primitive_c, Primitive_p;

    Conservative2Primitive(Primitive_mm, Conservative_mm, EOS_Parameter);
    Conservative2Primitive(Primitive_m, Conservative_m, EOS_Parameter);
    Conservative2Primitive(Primitive_c, Conservative_c, EOS_Parameter);
    Conservative2Primitive(Primitive_p, Conservative_p, EOS_Parameter);

    amrex::Array<amrex::Real, 4> Alpha = {
        Conservative_mm[INDEX_VolumeFraction1],
        Conservative_m[INDEX_VolumeFraction1],
        Conservative_c[INDEX_VolumeFraction1],
        Conservative_p[INDEX_VolumeFraction1]
    };
    

    bool Gradient_L = (Conservative_c[INDEX_VolumeFraction1] - Conservative_m[INDEX_VolumeFraction1])*(Conservative_m[INDEX_VolumeFraction1] - Conservative_mm[INDEX_VolumeFraction1]) > 0.0;
    bool Gradient_R = (Conservative_p[INDEX_VolumeFraction1] - Conservative_c[INDEX_VolumeFraction1])*(Conservative_c[INDEX_VolumeFraction1] - Conservative_m[INDEX_VolumeFraction1]) > 0.0;
    
    bool Gradient = Gradient_L && Gradient_R;
    bool Symmetric = (Conservative_mm[INDEX_VolumeFraction1] == Conservative_p[INDEX_VolumeFraction1]) && 
                     (Conservative_m[INDEX_VolumeFraction1] == Conservative_c[INDEX_VolumeFraction1]) &&
                     // Four equal values are not a two-cell plateau
                     (Conservative_m[INDEX_VolumeFraction1] != Conservative_mm[INDEX_VolumeFraction1]);

    // Without UseNormal, THINC is applied at every interface cell, with no normal check.
    bool NormalGate = (!ID_THINC_UseNormal) || (std::abs(Normal[ID_Dimension]) > NormalThreshold);

    bool Int_L = (Conservative_m[INDEX_VolumeFraction1] > EPS && Conservative_m[INDEX_VolumeFraction1] < 1.0 - EPS) && (Gradient_L || Symmetric) && NormalGate;
    bool Int_R = (Conservative_c[INDEX_VolumeFraction1] > EPS && Conservative_c[INDEX_VolumeFraction1] < 1.0 - EPS) && (Gradient_R || Symmetric) && NormalGate;

    if ((Int_L || Int_R)){
        FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                ConservativeL, ConservativeR,
                                                Conservative,
                                                EOS_Parameter,
                                                dx[ID_Dimension],
                                                ID_Dimension, 0, 1, 1);

        amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_L = ConservativeL;
        amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative_R = ConservativeR;


        amrex::Real beta = (ID_THINC_UseNormal) ? (Beta * std::abs(Normal[ID_Dimension]) + 0.01) : Beta;

        if (Int_L){
            amrex::Real sigma = SignNZ(Alpha[2] - Alpha[0]);
            amrex::Real Alpha1_Reconstructed = Reconstruction2Face_THINC_D(beta, sigma, Alpha[1], 1);
            Reconstruct_State_From_VolumeFraction(ConservativeL, Conservative_L, Alpha1_Reconstructed, EOS_Parameter);
        }

        if (Int_R){
            amrex::Real sigma = SignNZ(Alpha[3] - Alpha[1]);
            amrex::Real Alpha1_Reconstructed = Reconstruction2Face_THINC_D(beta, sigma, Alpha[2], 0);
            Reconstruct_State_From_VolumeFraction(ConservativeR, Conservative_R, Alpha1_Reconstructed, EOS_Parameter);        
        }

    } else {
        FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                ConservativeL, ConservativeR,
                                                Conservative,
                                                EOS_Parameter,
                                                dx[ID_Dimension],
                                                ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);
    }
}

#endif
#endif