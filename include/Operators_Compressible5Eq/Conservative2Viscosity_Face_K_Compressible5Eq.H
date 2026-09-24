// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Viscosity_Face_K_Compressible5Eq_H_
#define Conservative2Viscosity_Face_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Viscosity_Face_K(int i, int j, int k,
                                   amrex::Real& mu_Face,
                                   amrex::Real& mub_Face,
                                   amrex::Array4<amrex::Real const> const& Conservative,
                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                      AMREX_D_DECL(amrex::Real const& dx,
                                   amrex::Real const& dy,
                                   amrex::Real const& dz),
                                   amrex::Real const& time,
                                   Parm const& EOS_Parameter,
                                   int ID_Dimension)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_m;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_c;
    for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
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
    Conservative2Viscosity_Face_D(mu_Face,
                                  mub_Face,
                                  Conservative_m,
                                  Conservative_c,
                                  EOS_Parameter);
    
}

#endif
#endif