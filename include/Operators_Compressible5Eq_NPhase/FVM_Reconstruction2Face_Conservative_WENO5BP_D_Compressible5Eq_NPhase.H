// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO5BP_D_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_WENO5BP_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO5BP_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Conservative_Face,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_mm,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_m,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_p,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_pp,
                                                    Parm const& EOS_Parameter,
                                                    amrex::Real const& dx,
                                                    int ID_Face,
                                                    int ID_Poly)
{
    int ID_FaceRe;
    //fL_fp
    ID_FaceRe = 1;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeL_fp;
    FVM_Reconstruction2Face_Conservative_WENO5B_D(ConservativeL_fp,
                                                  Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp,
                                                  EOS_Parameter,
                                                  dx,
                                                  ID_FaceRe, ID_Poly);
    //fR_fm
    ID_FaceRe = -1;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeR_fm;
    FVM_Reconstruction2Face_Conservative_WENO5B_D(ConservativeR_fm,
                                                  Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp,
                                                  EOS_Parameter,
                                                  dx,
                                                  ID_FaceRe, ID_Poly);

    //EnergyBound
    amrex::Real Weight;
    if (ID_Poly == 2){
        Weight = 1.0/6.0;
    }
    else if (ID_Poly == 4){
        Weight = 1.0/12.0;
    }
    FVM_Reconstruction2Face_Conservative_LimiterBound(Conservative_Face,
                                                      Conservative_c, ConservativeL_fp, ConservativeR_fm,
                                                      EOS_Parameter,
                                                      Weight,
                                                      ID_Face);
    
}

#endif
#endif