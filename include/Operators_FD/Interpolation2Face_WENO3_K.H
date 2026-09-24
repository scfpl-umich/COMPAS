// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Face_WENO3_K_H_
#define Interpolation2Face_WENO3_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Face_WENO3_K(int i, int j, int k, int n,
                                amrex::Real& f_Face,
                                amrex::Array4<amrex::Real const> const& f,
                                amrex::Real const& V,
                                int ID_Dimension)
{
    amrex::Real f_mm;
    amrex::Real f_m;
    amrex::Real f_c;
    amrex::Real f_p;
    
    f_c = f(i,j,k,n);
    //x
    if (ID_Dimension == 0){
        f_mm = f(i-2,j,k,n);
        f_m = f(i-1,j,k,n);
        f_p = f(i+1,j,k,n);
    }
    //y
    else if (ID_Dimension == 1){
        f_mm = f(i,j-2,k,n);
        f_m = f(i,j-1,k,n);
        f_p = f(i,j+1,k,n);
    }
    //z
    else if (ID_Dimension == 2){
        f_mm = f(i,j,k-2,n);
        f_m = f(i,j,k-1,n);
        f_p = f(i,j,k+1,n);
    }
    
    Interpolation2Face_WENO3_D(f_Face,
                               f_mm,
                               f_m,
                               f_c,
                               f_p,
                               V);
    
}

#endif