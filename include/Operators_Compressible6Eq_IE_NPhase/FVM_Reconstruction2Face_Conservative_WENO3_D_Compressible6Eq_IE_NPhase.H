// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO3_D_Compressible6Eq_IE_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_WENO3_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO3_D(amrex::Array<amrex::Real,NSTATE>& Conservative_Face,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_m,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_c,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_p,
                                                  Parm const& EOS_Parameter,
                                                  amrex::Real const& dx,
                                                  int ID_Face)
{
    amrex::Array<amrex::Real,NSTATE> Primitive_m;
    Conservative2Primitive(Primitive_m,
                           Conservative_m,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_p;
    Conservative2Primitive(Primitive_p,
                           Conservative_p,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,NSTATE> Primitive_Face;
    for (int i_State = 0; i_State < NSTATE; i_State++){
        amrex::Real f_m = Primitive_m[i_State];
        amrex::Real f_c = Primitive_c[i_State];
        amrex::Real f_p = Primitive_p[i_State];
        
        amrex::Real f_Face;
        Reconstruction2Face_WENO3_D(f_Face,
                                    f_m, f_c, f_p,
                                    ID_Face);
        Primitive_Face[i_State] = f_Face;
    }

    Primitive2Conservative(Conservative_Face,
                        Primitive_Face,
                        EOS_Parameter);
    
    //Consistent limiter
    // amrex::Array<amrex::Real,NSTATE> Conservative_UnLimit;
    // Primitive2Conservative(Conservative_UnLimit,
    //                        Primitive_Face,
    //                        EOS_Parameter);
    // FVM_Reconstruction2Face_Conservative_LimiterConsistent(Conservative_Face,
    //                                                        Conservative_UnLimit,
    //                                                        Conservative_c,
    //                                                        EOS_Parameter,
    //                                                        dx,
    //                                                        ID_Face);
    
}

#endif
#endif