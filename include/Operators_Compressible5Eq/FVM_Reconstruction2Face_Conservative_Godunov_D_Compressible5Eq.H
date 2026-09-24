// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_Godunov_D_Compressible5Eq_H_
#define FVM_Reconstruction2Face_Conservative_Godunov_D_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_Godunov_D(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Conservative_Face,
                                                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_c,
                                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_Face;
    for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
        amrex::Real f_c  = Primitive_c [i_State];
        
        {
        BL_PROFILE("FVM_Reconstruction2Face_Conservative_D()::Reconstruction2Face_D");
        amrex::Real f_Face;
        Reconstruction2Face_Godunov_D(f_Face,
                                      f_c);
        Primitive_Face[i_State] = f_Face;
        }
    }
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif