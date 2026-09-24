// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_Godunov_D_Compressible6Eq_IE_NPhase_H_
#define FVM_Reconstruction2FaceLR_Conservative_Godunov_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_Godunov_D(amrex::Array<amrex::Real,NSTATE>& ConservativeL,
                                                      amrex::Array<amrex::Real,NSTATE>& ConservativeR,
                                                      amrex::Array<amrex::Real,NSTATE> const& Conservative_m,
                                                      amrex::Array<amrex::Real,NSTATE> const& Conservative_c,
                                                      Parm const& EOS_Parameter)
{
    int ID_Face;
    
    //fL_fm
    ID_Face = 1;
    FVM_Reconstruction2Face_Conservative_Godunov_D(ConservativeL,
                                                   Conservative_m,
                                                   EOS_Parameter);
    
    //fR_fm
    ID_Face = -1;
    FVM_Reconstruction2Face_Conservative_Godunov_D(ConservativeR,
                                                   Conservative_c,
                                                   EOS_Parameter);
    
}

#endif
#endif