// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_WENO5B_K_Compressible5Eq_H_
#define FVM_Reconstruction2FaceLR_Conservative_WENO5B_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_WENO5B_K(int i, int j, int k,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeL,
                                                     amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& ConservativeR,
                                                     amrex::Array4<amrex::Real const> const& Conservative,
                                                     Parm const& EOS_Parameter,
                                                     int ID_Dimension,
                                                     int ID_Poly)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_mmm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_mm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_m;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_c;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_p;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_pp;
    for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
        Conservative_c[i_State] = Conservative(i,j,k,i_State);
        //x
        if (ID_Dimension == 0){
            Conservative_mmm[i_State] = Conservative(i-3,j,k,i_State);
            Conservative_mm [i_State] = Conservative(i-2,j,k,i_State);
            Conservative_m  [i_State] = Conservative(i-1,j,k,i_State);
            Conservative_p  [i_State] = Conservative(i+1,j,k,i_State);
            Conservative_pp [i_State] = Conservative(i+2,j,k,i_State);
        }
        //y
        else if (ID_Dimension == 1){
            Conservative_mmm[i_State] = Conservative(i,j-3,k,i_State);
            Conservative_mm [i_State] = Conservative(i,j-2,k,i_State);
            Conservative_m  [i_State] = Conservative(i,j-1,k,i_State);
            Conservative_p  [i_State] = Conservative(i,j+1,k,i_State);
            Conservative_pp [i_State] = Conservative(i,j+2,k,i_State);
        }
        //z
        else if (ID_Dimension == 2){
            Conservative_mmm[i_State] = Conservative(i,j,k-3,i_State);
            Conservative_mm [i_State] = Conservative(i,j,k-2,i_State);
            Conservative_m  [i_State] = Conservative(i,j,k-1,i_State);
            Conservative_p  [i_State] = Conservative(i,j,k+1,i_State);
            Conservative_pp [i_State] = Conservative(i,j,k+2,i_State);
        }
    }
    FVM_Reconstruction2FaceLR_Conservative_WENO5B_D(ConservativeL, ConservativeR,
                                                    Conservative_mmm, Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp,
                                                    EOS_Parameter,
                                                    ID_Poly);
    
}

#endif
#endif