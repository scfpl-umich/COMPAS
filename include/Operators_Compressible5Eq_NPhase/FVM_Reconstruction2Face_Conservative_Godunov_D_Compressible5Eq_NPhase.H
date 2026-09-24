// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_Godunov_D_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_Godunov_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_Godunov_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Conservative_Face,
                                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_Face;
    for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
        amrex::Real f_c  = Primitive_c [i_State];
        
        amrex::Real f_Face;
        Reconstruction2Face_Godunov_D(f_Face,
                                      f_c);
        Primitive_Face[i_State] = f_Face;
    }
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif