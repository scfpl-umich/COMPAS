// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_THINC_K_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_THINC_K_Compressible5Eq_NPhase_H_

#include<Macros.H>


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
Reconstruct_State_From_VolumeFraction (amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> & Conservative_Reconstructed,
                                       amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const & Conservative,
                                       amrex::Array<amrex::Real,NPHASE> & Alpha_Reconstructed,
                                       Parm const& EOS_Parameter)
{

    amrex::Real eps = 1.0e-12;    

    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;)
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;

    Conservative2Primitive(Primitive,Conservative,EOS_Parameter);
    Primitive2State(Alpha_rho, AMREX_D_DECL(u,v,w), p, Alpha, Primitive);
    
    amrex::Real K = 0.5 * (AMREX_D_TERM(u*u, + v*v, + w*w));
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Recon;

    amrex::Real rho_Recon = 0.0;
    amrex::Real rho_E_Recon = 0.0;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){


        amrex::Real rho_k = Alpha_rho[i_Phase] / (Alpha[i_Phase] + 1.0e-16);
        amrex::Real rho_e_k = EOS_Parameter.A[i_Phase] * p + EOS_Parameter.B[i_Phase] + rho_k * EOS_Parameter.D[i_Phase];

        Alpha_rho_Recon[i_Phase] = Alpha_Reconstructed[i_Phase] * rho_k;

        rho_Recon += Alpha_rho_Recon[i_Phase];
        rho_E_Recon += Alpha_Reconstructed[i_Phase] * rho_e_k;
    }

    rho_E_Recon += rho_Recon * K;

    AMREX_D_TERM(
    amrex::Real Mom_x_Recon = u * rho_Recon;,
    amrex::Real Mom_y_Recon = v * rho_Recon;,
    amrex::Real Mom_z_Recon = w * rho_Recon;)
    
    Component2Vector(Conservative_Reconstructed, 
                     Alpha_rho_Recon,
                     AMREX_D_DECL(
                     u * rho_Recon,
                     v * rho_Recon,
                     w * rho_Recon),
                     rho_E_Recon,
                     Alpha_Reconstructed);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
FVM_Reconstruction2Face_Conservative_THINC_K(int i, int j, int k,
                                             amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeL,
                                             amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeR,
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

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;

    for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
        Conservative_c[i_State] = Conservative(i,j,k,i_State);
        //x
        if (ID_Dimension == 0){
            Conservative_mm[i_State] = Conservative(i-2,j,k,i_State);
            Conservative_m [i_State] = Conservative(i-1,j,k,i_State);
            Conservative_p [i_State] = Conservative(i+1,j,k,i_State);
        }
        //y
        else if (ID_Dimension == 1){
            Conservative_mm[i_State] = Conservative(i,j-2,k,i_State);
            Conservative_m [i_State] = Conservative(i,j-1,k,i_State);
            Conservative_p [i_State] = Conservative(i,j+1,k,i_State);
        }
        //z
        else if (ID_Dimension == 2){
            Conservative_mm[i_State] = Conservative(i,j,k-2,i_State);
            Conservative_m [i_State] = Conservative(i,j,k-1,i_State);
            Conservative_p [i_State] = Conservative(i,j,k+1,i_State);
        }
    }

    // Primitive conversion if needed
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_mm, Primitive_m, Primitive_c, Primitive_p;
    Conservative2Primitive(Primitive_mm, Conservative_mm, EOS_Parameter);
    Conservative2Primitive(Primitive_m, Conservative_m, EOS_Parameter);
    Conservative2Primitive(Primitive_c, Conservative_c, EOS_Parameter);
    Conservative2Primitive(Primitive_p, Conservative_p, EOS_Parameter);

    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>, NPHASE> Normal;
    amrex::Array<amrex::Real, NPHASE> beta;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        int vidx = INDEX_VolumeFraction1 + i_Phase;
        if (ID_Dimension == 0){
            FDM_Normal_XFaceFromCenter_K(i,j,k,vidx, Normal[i_Phase], Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));
        }
        else if (ID_Dimension == 1){
#if (AMREX_SPACEDIM > 1)
            FDM_Normal_YFaceFromCenter_K(i,j,k,vidx, Normal[i_Phase], Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
        else if (ID_Dimension == 2){
#if (AMREX_SPACEDIM > 2)
            FDM_Normal_ZFaceFromCenter_K(i,j,k,vidx, Normal[i_Phase], Conservative, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
        beta[i_Phase] = (ID_THINC_UseNormal) ? (Beta * std::abs(Normal[i_Phase][ID_Dimension]) + 0.01) : Beta;
    }

    // THINC eligibility flags for left, right interface, and global
    amrex::Array<bool, NPHASE> THINC_used_L;
    amrex::Array<bool, NPHASE> THINC_used_R;
    bool any_THINC_used_L = false, any_THINC_used_R = false;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++) {
        int vidx = INDEX_VolumeFraction1 + i_Phase;

        amrex::Real a_mm = Conservative_mm[vidx];
        amrex::Real a_m  = Conservative_m[vidx];
        amrex::Real a_c  = Conservative_c[vidx];
        amrex::Real a_p  = Conservative_p[vidx];
        bool Symmetric = (a_mm == a_p) && (a_m == a_c) && (a_m != a_mm); // Four equal values are not a two-cell plateau
        // Without UseNormal, THINC is applied at every interface cell, with no normal check.
        bool NormalGate = (!ID_THINC_UseNormal) || (std::abs(Normal[i_Phase][ID_Dimension]) > NormalThreshold);

        // LEFT interface gates
        bool Interface_L = (a_m > EPS && a_m < 1.0 - EPS);
        bool Gradient_L = (a_c - a_m)*(a_m - a_mm) > 0.0;
        bool Int_L = Interface_L && (Gradient_L || Symmetric) && NormalGate;
        THINC_used_L[i_Phase] = Int_L;
        any_THINC_used_L |= Int_L;

        // RIGHT interface gates
        bool Interface_R = (a_c > EPS && a_c < 1.0 - EPS);
        bool Gradient_R = (a_p - a_c)*(a_c - a_m) > 0.0;
        bool Int_R = Interface_R && (Gradient_R || Symmetric) && NormalGate;
        THINC_used_R[i_Phase] = Int_R;
        any_THINC_used_R |= Int_R;
    }

    // THINC global
    bool use_THINC = (any_THINC_used_L || any_THINC_used_R);

    // --- Fallback: if no phase qualifies for THINC anywhere, use statically defined method ---
    if (use_THINC) {
        FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                ConservativeL, ConservativeR,
                                                Conservative,
                                                EOS_Parameter,
                                                dx[ID_Dimension],
                                                ID_Dimension, 0, 1, 1); 
    } else {
        FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                ConservativeL, ConservativeR,
                                                Conservative,
                                                EOS_Parameter,
                                                dx[ID_Dimension],
                                                ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);
        return;
    }

    // --- If THINC applies anywhere, reconstruct per phase and interface ---
    amrex::Array<amrex::Real,NPHASE> Alpha_Reconstructed;

    // LEFT interface
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_L = ConservativeL; // fallback starting point
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++) {
        int vidx = INDEX_VolumeFraction1 + i_Phase;
        amrex::Real a_mm = Conservative_mm[vidx];
        amrex::Real a_m  = Conservative_m[vidx];
        amrex::Real a_c  = Conservative_c[vidx];
        amrex::Real a_p  = Conservative_p[vidx];
        if (THINC_used_L[i_Phase]) {
            amrex::Real sigma = SignNZ(a_c - a_mm);
            Alpha_Reconstructed[i_Phase] = Reconstruction2Face_THINC_D(beta[i_Phase], sigma, a_m, 1);
        } else {
            Alpha_Reconstructed[i_Phase] = a_m; // fallback is just a_m
        }
    }
    Reconstruct_State_From_VolumeFraction(ConservativeL, Conservative_L, Alpha_Reconstructed, EOS_Parameter);  
    Conservative_L = ConservativeL;
    FVM_Reconstruction2Face_Conservative_LimiterConsistent(ConservativeL, Conservative_L, Conservative_m, EOS_Parameter, AMREX_D_PICK(dx[0],dx[1],dx[2]), 1);

    // RIGHT interface
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++) {
        int vidx = INDEX_VolumeFraction1 + i_Phase;
        amrex::Real a_mm = Conservative_mm[vidx];
        amrex::Real a_m  = Conservative_m[vidx];
        amrex::Real a_c  = Conservative_c[vidx];
        amrex::Real a_p  = Conservative_p[vidx];
        if (THINC_used_R[i_Phase]) {
            amrex::Real sigma = SignNZ(a_p - a_m);
            Alpha_Reconstructed[i_Phase] = Reconstruction2Face_THINC_D(beta[i_Phase], sigma, a_c, 0);
        } else {
            Alpha_Reconstructed[i_Phase] = a_c; // fallback is just a_c
        }
    }
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_R = ConservativeR;
    Reconstruct_State_From_VolumeFraction(ConservativeR, Conservative_R, Alpha_Reconstructed, EOS_Parameter);
    Conservative_R = ConservativeR;
    FVM_Reconstruction2Face_Conservative_LimiterConsistent(ConservativeR, Conservative_R, Conservative_c, EOS_Parameter, AMREX_D_PICK(dx[0],dx[1],dx[2]), -1);
}

#endif