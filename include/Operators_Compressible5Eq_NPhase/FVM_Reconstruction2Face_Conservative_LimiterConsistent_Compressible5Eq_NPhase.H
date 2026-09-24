// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_LimiterConsistent_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_LimiterConsistent_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_LimiterConsistent(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Conservative_Limit,
                                                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_Face,
                                                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                                            Parm const& EOS_Parameter,
                                                            amrex::Real const& dx,
                                                            int ID_Face)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_Face;
    Conservative2Primitive(Primitive_Face,
                           Conservative_Face,
                           EOS_Parameter);

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Face;
    AMREX_D_TERM(amrex::Real u_Face;,
                 amrex::Real v_Face;,
                 amrex::Real w_Face);
                 amrex::Real p_Face;
    amrex::Array<amrex::Real,NPHASE> Alpha_Face;
    Primitive2State(Alpha_rho_Face,
       AMREX_D_DECL(u_Face, v_Face, w_Face),
                    p_Face,
                    Alpha_Face,
                    Primitive_Face);

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_c;
    AMREX_D_TERM(amrex::Real u_c;,
                 amrex::Real v_c;,
                 amrex::Real w_c);
                 amrex::Real p_c;
    amrex::Array<amrex::Real,NPHASE> Alpha_c;
    Primitive2State(Alpha_rho_c,
       AMREX_D_DECL(u_c, v_c, w_c),
                    p_c,
                    Alpha_c,
                    Primitive_c);

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Limit;
    amrex::Array<amrex::Real,NPHASE> Alpha_Limit;
    FVM_Reconstruction2Face_LimiterConsistent(Alpha_rho_Limit, Alpha_Limit,
                                              Alpha_rho_Face , Alpha_Face,
                                              Alpha_rho_c    , Alpha_c,
                                              dx,
                                              ID_Face);

    State2Primitive(Primitive_Face,
                    Alpha_rho_Limit,
       AMREX_D_DECL(u_Face, v_Face, w_Face),
                    p_Face,
                    Alpha_Limit);
    Primitive2Conservative(Conservative_Limit,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif