// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_THINCOptions_K_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_THINCOptions_K_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

// N-phase THINC with FiniteVolume.THINC.Multi, FaceState and PhaseDensity, used when one of them is
// not 0; otherwise FVM_Reconstruction2Face_Conservative_THINC_K. Each phase is selected on its own.
//   FaceState = 0:    both sides of a face with a flagged phase are rebuilt from MUSCL-minmod; the
//                     phases that are not flagged take their cell average.
//   FaceState = 1:    the base scheme for every variable. Only a side with a flagged phase is
//                     rebuilt: flagged phases take their THINC value, the others keep their face
//                     value, and the face volume fractions are divided by their sum (Zhang et al.
//                     2025, Eq. 47).
//   PhaseDensity = 1: rho_k of a rebuilt side is the phase density of the cell,
//                     (alpha_k rho_k)_cell / alpha_k,cell (Shyue and Xiao 2014), instead of the
//                     ratio of the face values.
// The consistent limiter is applied to every rebuilt side.

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
Reconstruct_State_From_VolumeFraction_Options (amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> & Conservative_Reconstructed,
                                               amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const & Conservative_Face,
                                               amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const & Conservative_Cell,
                                               amrex::Array<amrex::Real,NPHASE> const & Alpha_Reconstructed,
                                               int ID_PhaseDensity,
                                               int ID_FaceState,
                                               Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;)
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;

    Conservative2Primitive(Primitive,Conservative_Face,EOS_Parameter);
    Primitive2State(Alpha_rho, AMREX_D_DECL(u,v,w), p, Alpha, Primitive);

    amrex::Real K = 0.5 * (AMREX_D_TERM(u*u, + v*v, + w*w));
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Recon;
    amrex::Array<amrex::Real,NPHASE> Alpha_Recon;

    amrex::Real rho_Recon = 0.0;
    amrex::Real rho_E_Recon = 0.0;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real rho_k;
        // With FaceState = 1 the face value of alpha_k can be much smaller than the cell value,
        // and the ratio of the face values is then ill-conditioned: use the cell value there
        bool Cell_Density = (ID_PhaseDensity == 1) ||
                            (ID_FaceState == 1 && !(Alpha[i_Phase] >= 0.25 * Conservative_Cell[INDEX_VolumeFraction1 + i_Phase]));
        if (Cell_Density){
            rho_k = Conservative_Cell[INDEX_Mass1 + i_Phase] / (Conservative_Cell[INDEX_VolumeFraction1 + i_Phase] + 1.0e-16);
        } else {
            rho_k = Alpha_rho[i_Phase] / (Alpha[i_Phase] + 1.0e-16);
        }
        amrex::Real rho_e_k = EOS_Parameter.A[i_Phase] * p + EOS_Parameter.B[i_Phase] + rho_k * EOS_Parameter.D[i_Phase];
        // NASG: rho_k e_k = (1 - b_k rho_k)(A_k p + B_k) + rho_k q_k
        if (EOS_Parameter.C[i_Phase] > 0.0){
            rho_e_k = NASG_Xi_rho(rho_k, EOS_Parameter.C[i_Phase]) * (EOS_Parameter.A[i_Phase] * p + EOS_Parameter.B[i_Phase]) + rho_k * EOS_Parameter.D[i_Phase];
        }

        Alpha_Recon[i_Phase] = Alpha_Reconstructed[i_Phase];
        Alpha_rho_Recon[i_Phase] = Alpha_Reconstructed[i_Phase] * rho_k;

        rho_Recon += Alpha_rho_Recon[i_Phase];
        rho_E_Recon += Alpha_Reconstructed[i_Phase] * rho_e_k;
    }

    rho_E_Recon += rho_Recon * K;

    Component2Vector(Conservative_Reconstructed,
                     Alpha_rho_Recon,
                     AMREX_D_DECL(
                     u * rho_Recon,
                     v * rho_Recon,
                     w * rho_Recon),
                     rho_E_Recon,
                     Alpha_Recon);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
FVM_Reconstruction2Face_Conservative_THINCOptions_K(int i, int j, int k,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeL,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeR,
                                                    amrex::Array4<amrex::Real const> const& Conservative,
                                                    Parm const& EOS_Parameter,
                                                    amrex::Array<amrex::Real, AMREX_SPACEDIM> const& dx,
                                                    int ID_Dimension,
                                                    FiniteVolume_Parameter const& FV)
{
    constexpr int NS = 2*NPHASE+AMREX_SPACEDIM+1;

    amrex::Array<bool,NPHASE> Int_L, Int_R;
    amrex::Array<amrex::Real,NPHASE> Alpha_L, Alpha_R;
    bool any_L = false, any_R = false;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        bool iL, iR;
        THINC_Select_K(i, j, k, ID_Dimension, INDEX_VolumeFraction1 + i_Phase, Conservative, dx, FV,
                       iL, iR, Alpha_L[i_Phase], Alpha_R[i_Phase]);
        Int_L[i_Phase] = iL;
        Int_R[i_Phase] = iR;
        any_L = any_L || iL;
        any_R = any_R || iR;
    }

    bool Minmod = (any_L || any_R) && (FV.THINC_FaceState == 0);
    FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                            ConservativeL, ConservativeR,
                                            Conservative,
                                            EOS_Parameter,
                                            dx[ID_Dimension],
                                            ID_Dimension,
                                            (Minmod) ? 0 : FV.ID_WENO,
                                            (Minmod) ? 1 : FV.ID_Limiter,
                                            (Minmod) ? 1 : FV.ID_Poly);
    if (!any_L && !any_R){ return; }

    const int ei = (ID_Dimension == 0), ej = (ID_Dimension == 1), ek = (ID_Dimension == 2);
    amrex::Array<amrex::Real,NS> Conservative_m, Conservative_c;
    for (int n = 0; n < NS; n++){
        Conservative_m[n] = Conservative(i-ei, j-ej, k-ek, n);
        Conservative_c[n] = Conservative(i, j, k, n);
    }

    for (int Side = 0; Side < 2; Side++){
        bool any = (Side == 0) ? any_L : any_R;
        if (FV.THINC_FaceState == 1 && !any){ continue; }
        amrex::Array<amrex::Real,NS>& Face = (Side == 0) ? ConservativeL : ConservativeR;
        amrex::Array<amrex::Real,NS> const& Cell = (Side == 0) ? Conservative_m : Conservative_c;

        amrex::Array<amrex::Real,NPHASE> Alpha_Reconstructed;
        amrex::Real Sum = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            int vidx = INDEX_VolumeFraction1 + i_Phase;
            bool Flag = (Side == 0) ? Int_L[i_Phase] : Int_R[i_Phase];
            if (Flag){
                Alpha_Reconstructed[i_Phase] = (Side == 0) ? Alpha_L[i_Phase] : Alpha_R[i_Phase];
            } else {
                Alpha_Reconstructed[i_Phase] = (FV.THINC_FaceState == 0) ? Cell[vidx] : Face[vidx];
            }
            Sum += Alpha_Reconstructed[i_Phase];
        }
        if (FV.THINC_FaceState == 1){
            for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){ Alpha_Reconstructed[i_Phase] /= Sum; }
        }
        amrex::Array<amrex::Real,NS> Conservative_Face = Face;
        Reconstruct_State_From_VolumeFraction_Options(Face, Conservative_Face, Cell, Alpha_Reconstructed,
                                                      FV.THINC_PhaseDensity, FV.THINC_FaceState, EOS_Parameter);
        Conservative_Face = Face;
        FVM_Reconstruction2Face_Conservative_LimiterConsistent(Face, Conservative_Face, Cell, EOS_Parameter,
                                                               AMREX_D_PICK(dx[0],dx[1],dx[2]), (Side == 0) ? 1 : -1);
    }
}

#endif
#endif
