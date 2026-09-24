// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_WENO5_K_H_
#define Reconstruction_WENO5_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_WENO5_K(int i, int j, int k, int n,
                            amrex::Real& f_Reconstruction,
                            amrex::Real const& x,
                            amrex::Array4<amrex::Real const> const& f,
                            amrex::Real const& dx,
                            amrex::Real const& prob_lo,
                            int ID_Dimension,
                            int ID_Poly)
{
    int i_Index;
    amrex::Real f_mm;
    amrex::Real f_m;
    amrex::Real f_c;
    amrex::Real f_p;
    amrex::Real f_pp;
    
    f_c = f(i,j,k,n);
    //x
    if (ID_Dimension == 0){
        f_mm = f(i-2,j,k,n);
        f_m = f(i-1,j,k,n);
        f_p = f(i+1,j,k,n);
        f_pp = f(i+2,j,k,n);
        i_Index = i;
    }
    //y
    else if (ID_Dimension == 1){
        f_mm = f(i,j-2,k,n);
        f_m = f(i,j-1,k,n);
        f_p = f(i,j+1,k,n);
        f_pp = f(i,j+2,k,n);
        i_Index = j;
    }
    //z
    else if (ID_Dimension == 2){
        f_mm = f(i,j,k-2,n);
        f_m = f(i,j,k-1,n);
        f_p = f(i,j,k+1,n);
        f_pp = f(i,j,k+2,n);
        i_Index = k;
    }
	amrex::Real x_c = prob_lo + ( 0.5 + amrex::Real(i_Index) ) * dx;
    
    Reconstruction_WENO5_D(f_Reconstruction,
                           x,
                           f_mm,
                           f_m,
                           f_c,
                           f_p,
                           f_pp,
                           x_c,
                           dx,
                           ID_Poly);
    
}

#endif