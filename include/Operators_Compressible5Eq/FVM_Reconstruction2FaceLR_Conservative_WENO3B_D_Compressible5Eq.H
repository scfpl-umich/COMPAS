// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_WENO3B_D_Compressible5Eq_H_
#define FVM_Reconstruction2FaceLR_Conservative_WENO3B_D_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_WENO3B_D(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeL,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeR,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_mm,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_m,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_c,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_p,
                                                     Parm const& EOS_Parameter)
{
    int ID_Face;
    {
    BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::FVM_Reconstruction2Face_Conservative_D");
    //fL_fm
    ID_Face = 1;
    FVM_Reconstruction2Face_Conservative_WENO3B_D(ConservativeL,
                                                  Conservative_mm, Conservative_m, Conservative_c,
                                                  EOS_Parameter,
                                                  ID_Face);
    
    //fR_fm
    ID_Face = -1;
    FVM_Reconstruction2Face_Conservative_WENO3B_D(ConservativeR,
                                                  Conservative_m, Conservative_c, Conservative_p,
                                                  EOS_Parameter,
                                                  ID_Face);
    }
    
}

#endif
#endif