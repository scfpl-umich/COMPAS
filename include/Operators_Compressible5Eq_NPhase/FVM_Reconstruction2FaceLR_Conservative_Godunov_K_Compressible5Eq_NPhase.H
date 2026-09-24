// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_Godunov_K_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2FaceLR_Conservative_Godunov_K_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_Godunov_K(int i, int j, int k,
                                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeL,
                                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeR,
                                                      amrex::Array4<amrex::Real const> const& Conservative,
                                                      Parm const& EOS_Parameter,
                                                      int ID_Dimension)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
    for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
        Conservative_c[i_State] = Conservative(i,j,k,i_State);
        //x
        if (ID_Dimension == 0){
            Conservative_m[i_State] = Conservative(i-1,j,k,i_State);
        }
        //y
        else if (ID_Dimension == 1){
            Conservative_m[i_State] = Conservative(i,j-1,k,i_State);
        }
        //z
        else if (ID_Dimension == 2){
            Conservative_m[i_State] = Conservative(i,j,k-1,i_State);
        }
    }
    FVM_Reconstruction2FaceLR_Conservative_Godunov_D(ConservativeL, ConservativeR,
                                                     Conservative_m, Conservative_c,
                                                     EOS_Parameter);
    
}

#endif
#endif