// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_CDI_K_H_
#define FDM_PhaseFieldFlux_CDI_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_CDI_K(int i, int j, int k, int n,
                              amrex::Real& PF_Flux,
                              amrex::Array4<amrex::Real const> const& PF,
                              amrex::Array4<amrex::Real const> const& Normal,
                              amrex::Array4<amrex::Real const> const& Mob,
                              amrex::Real const& Eta,
                              amrex::Real const& dx,
                              int ID_Dimension)
{
    amrex::Real PF_c = PF(i,j,k,n);
    amrex::Real Normal_c = Normal(i,j,k,ID_Dimension);
    
    amrex::Real PF_m;
    amrex::Real Normal_m;
    //x
    if (ID_Dimension == 0){
        PF_m = PF(i-1,j,k,n);
        Normal_m = Normal(i-1,j,k,ID_Dimension);
    }
    //y
    else if (ID_Dimension == 1){
        PF_m = PF(i,j-1,k,n);
        Normal_m = Normal(i,j-1,k,ID_Dimension);
    }
    //z
    else if (ID_Dimension == 2){
        PF_m = PF(i,j,k-1,n);
        Normal_m = Normal(i,j,k-1,ID_Dimension);
    }
    
    amrex::Real Mob_Face;
    Interpolation2Face_Linear_K(i, j, k, 0,
                                Mob_Face,
                                Mob,
                                ID_Dimension);
    
    FDM_PhaseFieldFlux_CDI_D(PF_Flux,
                             PF_c, PF_m,
                             Normal_c, Normal_m,
                             Mob_Face, Eta,
                             dx);
}

#endif