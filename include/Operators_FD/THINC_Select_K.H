// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef THINC_Select_K_H_
#define THINC_Select_K_H_

// THINC sides and THINC face values of a volume fraction for FiniteVolume.THINC.Multi, FaceState
// and PhaseDensity (FVM_Reconstruction2Face_Conservative_THINCOptions_K). The selection is that of
// FVM_Reconstruction2Face_Conservative_THINC_K.

// Unit normal of cell (ii,jj,kk) in index space, pointing to larger values of component comp:
// Youngs' gradient (3x3 or 3x3x3 block, weights 1 2 1) of the distance transform of the volume
// fraction (Paula et al. 2023). Returns false where the gradient vanishes.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool THINC_CellNormal_K(int ii, int jj, int kk, int comp,
                        amrex::Array4<amrex::Real const> const& f,
                        amrex::Real Beta,
                        amrex::Array<amrex::Real,AMREX_SPACEDIM>& Normal)
{
    const amrex::Real d_corr = std::log(std::sinh(Beta) / Beta) / (2.0 * Beta);
    auto F = [=] (int a, int b, int c) -> amrex::Real
    {
        return THINC_DistanceTransform_D(f(ii+a, jj+b, kk+c, comp), Beta, d_corr);
    };
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Grad;
    for (int d = 0; d < AMREX_SPACEDIM; d++){ Grad[d] = 0.0; }
    const amrex::Real w[3] = {1.0, 2.0, 1.0};
#if (AMREX_SPACEDIM == 1)
    Grad[0] = 0.5 * (F(1,0,0) - F(-1,0,0));
#elif (AMREX_SPACEDIM == 2)
    for (int b = -1; b <= 1; b++){
        for (int a = -1; a <= 1; a++){
            if (a == 0 && b == 0){ continue; }
            amrex::Real v = F(a,b,0);
            Grad[0] += amrex::Real(a) * w[b+1] * v;
            Grad[1] += amrex::Real(b) * w[a+1] * v;
        }
    }
    Grad[0] *= 0.125;
    Grad[1] *= 0.125;
#else
    for (int c = -1; c <= 1; c++){
        for (int b = -1; b <= 1; b++){
            for (int a = -1; a <= 1; a++){
                if (a == 0 && b == 0 && c == 0){ continue; }
                amrex::Real v = F(a,b,c);
                Grad[0] += amrex::Real(a) * w[b+1] * w[c+1] * v;
                Grad[1] += amrex::Real(b) * w[a+1] * w[c+1] * v;
                Grad[2] += amrex::Real(c) * w[a+1] * w[b+1] * v;
            }
        }
    }
    Grad[0] /= 32.0;
    Grad[1] /= 32.0;
    Grad[2] /= 32.0;
#endif
    amrex::Real Sum = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; d++){ Sum += Grad[d] * Grad[d]; }
    if (Sum == 0.0){ return false; }
    amrex::Real Inv = 1.0 / std::sqrt(Sum);
    for (int d = 0; d < AMREX_SPACEDIM; d++){ Normal[d] = Grad[d] * Inv; }
    return true;
}

// THINC face value of cell (ii,jj,kk) at its right face (ID_Face = 1) or left face (ID_Face = 0);
// a_m, a_c, a_p are the cell and its neighbours along the direction.
// Multi = 0: 1D THINC (Reconstruction2Face_THINC_D) with the given beta.
// Multi = 1: face average of the MTHINC profile of the cell, with the constant Beta.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real THINC_FaceValue_K(int ii, int jj, int kk, int comp, int ID_Dimension, int ID_Face,
                              amrex::Real beta,
                              amrex::Real a_m, amrex::Real a_c, amrex::Real a_p,
                              amrex::Array4<amrex::Real const> const& f,
                              FiniteVolume_Parameter const& FV)
{
    if (FV.THINC_Multi == 0){
        return Reconstruction2Face_THINC_D(beta, SignNZ(a_p - a_m), a_c, ID_Face);
    }
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal;
    if (!THINC_CellNormal_K(ii, jj, kk, comp, f, FV.THINC_Beta, Normal)){
        return a_c;
    }
    amrex::Real d = MTHINC_Solve_D(FV.THINC_Beta, Normal, a_c);
    return MTHINC_FaceAverage_D(FV.THINC_Beta, Normal, d, ID_Dimension, (ID_Face == 1) ? 1 : -1);
}

// THINC flags and face values of component comp at the face between cell (i,j,k) - e_dim (left
// side) and cell (i,j,k) (right side): EPS window, monotone or Symmetric, and the |n_d| gate on
// the face normal.
//   Int_L, Int_R      THINC is used on the left / right side of the face
//   Alpha_L, Alpha_R  the THINC face value of that side (set only when the side is flagged)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void THINC_Select_K(int i, int j, int k, int ID_Dimension, int comp,
                    amrex::Array4<amrex::Real const> const& f,
                    amrex::Array<amrex::Real,AMREX_SPACEDIM> const& dx,
                    FiniteVolume_Parameter const& FV,
                    bool& Int_L, bool& Int_R,
                    amrex::Real& Alpha_L, amrex::Real& Alpha_R)
{
    const int ei = (ID_Dimension == 0), ej = (ID_Dimension == 1), ek = (ID_Dimension == 2);
    const amrex::Real EPS = FV.THINC_EPS;

    Int_L = false;
    Int_R = false;
    Alpha_L = 0.0;
    Alpha_R = 0.0;

    const amrex::Real a_mm = f(i-2*ei, j-2*ej, k-2*ek, comp);
    const amrex::Real a_m  = f(i-ei,   j-ej,   k-ek,   comp);
    const amrex::Real a_c  = f(i,      j,      k,      comp);
    const amrex::Real a_p  = f(i+ei,   j+ej,   k+ek,   comp);

    bool Window_L = (a_m > EPS && a_m < 1.0 - EPS);
    bool Window_R = (a_c > EPS && a_c < 1.0 - EPS);
    if (!Window_L && !Window_R){ return; }

    bool Gradient_L = (a_c - a_m) * (a_m - a_mm) > 0.0;
    bool Gradient_R = (a_p - a_c) * (a_c - a_m) > 0.0;
    // Four equal values are not a two-cell plateau
    bool Symmetric = (a_mm == a_p) && (a_m == a_c) && (a_m != a_mm);

    // Face normal (FDM_Normal_{X,Y,Z}FaceFromCenter_K) for the gate and for beta of the 1D THINC
    bool NormalGate = true;
    amrex::Real beta = FV.THINC_Beta;
    if (FV.ID_THINC_UseNormal){
        amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal;
        if (ID_Dimension == 0){
            FDM_Normal_XFaceFromCenter_K(i,j,k,comp, Normal, f, AMREX_D_DECL(dx[0],dx[1],dx[2]));
        }
        else if (ID_Dimension == 1){
#if (AMREX_SPACEDIM > 1)
            FDM_Normal_YFaceFromCenter_K(i,j,k,comp, Normal, f, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
        else if (ID_Dimension == 2){
#if (AMREX_SPACEDIM > 2)
            FDM_Normal_ZFaceFromCenter_K(i,j,k,comp, Normal, f, AMREX_D_DECL(dx[0],dx[1],dx[2]));
#endif
        }
        NormalGate = (std::abs(Normal[ID_Dimension]) > FV.THINC_NormalThreshold);
        beta = FV.THINC_Beta * std::abs(Normal[ID_Dimension]) + 0.01;
    }

    Int_L = Window_L && (Gradient_L || Symmetric) && NormalGate;
    Int_R = Window_R && (Gradient_R || Symmetric) && NormalGate;
    if (Int_L){
        Alpha_L = THINC_FaceValue_K(i-ei, j-ej, k-ek, comp, ID_Dimension, 1, beta, a_mm, a_m, a_c, f, FV);
    }
    if (Int_R){
        Alpha_R = THINC_FaceValue_K(i, j, k, comp, ID_Dimension, 0, beta, a_m, a_c, a_p, f, FV);
    }
}

#endif
