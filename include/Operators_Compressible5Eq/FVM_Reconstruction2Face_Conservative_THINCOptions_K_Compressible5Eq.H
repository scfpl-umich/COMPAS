// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_THINCOptions_K_Compressible5Eq_H_
#define FVM_Reconstruction2Face_Conservative_THINCOptions_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

// THINC with FiniteVolume.THINC.Multi, FaceState and PhaseDensity, used when one of them is not 0
// (FiniteVolume_Parameter::ID_THINC_Options = 1); otherwise FVM_Reconstruction2Face_Conservative_THINC_K.
//
// Face state of a THINC side from its THINC volume fraction Alpha1_Face and the face state
// Conservative_Face of the other variables (MUSCL-minmod, or the base scheme with FaceState = 1):
//   u, p from Conservative_Face;
//   rho_k = (alpha_k rho_k)_face / alpha_k,face    (PhaseDensity = 0; the cell value is used where
//           alpha_k,face < alpha_k,cell / 4, which the face values of the base scheme allow)
//   rho_k = (alpha_k rho_k)_cell / alpha_k,cell    (PhaseDensity = 1, Shyue and Xiao 2014)
//   (alpha_k rho_k) = alpha_k,THINC rho_k, E from pressure equilibrium.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
Reconstruct_State_From_VolumeFraction_Options (amrex::Array<amrex::Real,4+AMREX_SPACEDIM> & Conservative_Reconstructed,
                                               amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const & Conservative_Face,
                                               amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const & Conservative_Cell,
                                               amrex::Real Alpha1_Face,
                                               int ID_PhaseDensity,
                                               Parm const& EOS_Parameter)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;)
    amrex::Real p;
    amrex::Real Alpha1;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive;

    Conservative2Primitive(Primitive,Conservative_Face,EOS_Parameter);
    Primitive2State(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u,v,w), p, Alpha1, Primitive);

    amrex::Real rho1, rho2;
    if (ID_PhaseDensity == 1){
        amrex::Real Alpha1_Cell = Conservative_Cell[INDEX_VolumeFraction1];
        rho1 = Conservative_Cell[INDEX_Mass1] / Alpha1_Cell;
        rho2 = Conservative_Cell[INDEX_Mass2] / (1.0 - Alpha1_Cell);
    } else {
        rho1 = Alpha1_rho1 / Alpha1;
        rho2 = Alpha2_rho2 / (1.0 - Alpha1);
        amrex::Real Alpha1_Cell = Conservative_Cell[INDEX_VolumeFraction1];
        if (!(Alpha1 >= 0.25 * Alpha1_Cell)){             rho1 = Conservative_Cell[INDEX_Mass1] / Alpha1_Cell; }
        if (!(1.0 - Alpha1 >= 0.25 * (1.0 - Alpha1_Cell))){ rho2 = Conservative_Cell[INDEX_Mass2] / (1.0 - Alpha1_Cell); }
    }

    amrex::Real Alpha1_rho1_Recon = Alpha1_Face * rho1;
    amrex::Real Alpha2_rho2_Recon = (1.0 - Alpha1_Face) * rho2;
    amrex::Real rho_Recon = Alpha1_rho1_Recon + Alpha2_rho2_Recon;

    amrex::Real K = 0.5 * (AMREX_D_TERM(u*u, + v*v, + w*w));
    amrex::Real rho1_e1 = EOS_Parameter.A1 * p + EOS_Parameter.B1 + rho1 * EOS_Parameter.D1;
    amrex::Real rho2_e2 = EOS_Parameter.A2 * p + EOS_Parameter.B2 + rho2 * EOS_Parameter.D2;
    // NASG: rho_k e_k = (1 - b_k rho_k)(A_k p + B_k) + rho_k q_k
    if (EOS_Parameter.C1 > 0.0){ rho1_e1 = NASG_Xi_rho(rho1, EOS_Parameter.C1) * (EOS_Parameter.A1 * p + EOS_Parameter.B1) + rho1 * EOS_Parameter.D1; }
    if (EOS_Parameter.C2 > 0.0){ rho2_e2 = NASG_Xi_rho(rho2, EOS_Parameter.C2) * (EOS_Parameter.A2 * p + EOS_Parameter.B2) + rho2 * EOS_Parameter.D2; }

    Conservative_Reconstructed[INDEX_Mass1] = Alpha1_rho1_Recon;
    Conservative_Reconstructed[INDEX_Mass2] = Alpha2_rho2_Recon;
    AMREX_D_TERM(
    Conservative_Reconstructed[INDEX_MomentumX] = u * rho_Recon;,
    Conservative_Reconstructed[INDEX_MomentumY] = v * rho_Recon;,
    Conservative_Reconstructed[INDEX_MomentumZ] = w * rho_Recon;)
    Conservative_Reconstructed[INDEX_Energy] = K * rho_Recon + rho1_e1 * Alpha1_Face + rho2_e2 * (1.0 - Alpha1_Face);
    Conservative_Reconstructed[INDEX_VolumeFraction1] = Alpha1_Face;
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
FVM_Reconstruction2Face_Conservative_THINCOptions_K(int i, int j, int k,
                                                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeL,
                                                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeR,
                                                    amrex::Array4<amrex::Real const> const& Conservative,
                                                    Parm const& EOS_Parameter,
                                                    amrex::Array<amrex::Real, AMREX_SPACEDIM> const& dx,
                                                    int ID_Dimension,
                                                    FiniteVolume_Parameter const& FV)
{
    bool Int_L, Int_R;
    amrex::Real Alpha1_L, Alpha1_R;
    THINC_Select_K(i, j, k, ID_Dimension, INDEX_VolumeFraction1, Conservative, dx, FV,
                   Int_L, Int_R, Alpha1_L, Alpha1_R);

    // FaceState = 0: MUSCL-minmod for every variable on both sides of a THINC face
    // FaceState = 1: the base scheme for every variable; only the THINC side is rebuilt
    bool Minmod = (Int_L || Int_R) && (FV.THINC_FaceState == 0);
    FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                            ConservativeL, ConservativeR,
                                            Conservative,
                                            EOS_Parameter,
                                            dx[ID_Dimension],
                                            ID_Dimension,
                                            (Minmod) ? 0 : FV.ID_WENO,
                                            (Minmod) ? 1 : FV.ID_Limiter,
                                            (Minmod) ? 1 : FV.ID_Poly);
    if (!Int_L && !Int_R){ return; }

    const int ei = (ID_Dimension == 0), ej = (ID_Dimension == 1), ek = (ID_Dimension == 2);
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_Cell;
    if (Int_L){
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_Face = ConservativeL;
        for (int n = 0; n < 4+AMREX_SPACEDIM; n++){ Conservative_Cell[n] = Conservative(i-ei, j-ej, k-ek, n); }
        Reconstruct_State_From_VolumeFraction_Options(ConservativeL, Conservative_Face, Conservative_Cell,
                                                      Alpha1_L, FV.THINC_PhaseDensity, EOS_Parameter);
    }
    if (Int_R){
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_Face = ConservativeR;
        for (int n = 0; n < 4+AMREX_SPACEDIM; n++){ Conservative_Cell[n] = Conservative(i, j, k, n); }
        Reconstruct_State_From_VolumeFraction_Options(ConservativeR, Conservative_Face, Conservative_Cell,
                                                      Alpha1_R, FV.THINC_PhaseDensity, EOS_Parameter);
    }
}

#endif
#endif
