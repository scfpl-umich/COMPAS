// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Face_K_H_
#define Interpolation2Face_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Face_K(int i, int j, int k, int n,
                          amrex::Real& f_Face,
                          amrex::Array4<amrex::Real const> const& f,
                          amrex::Real const& f_max,
                          amrex::Real const& f_min,
                          amrex::Real const& V,
                          int ID_Dimension,
                          int ID_WENO,
                          int ID_Limiter,
                          int ID_Poly,
                          int ID_BoundType)
{
    //non-WENO
    if (ID_WENO == 0){
        //Linear (2nd order)
        if (ID_Limiter == 2){
            Interpolation2Face_Linear_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        ID_Dimension);
        }
        //Quick (3rd order)
        else if (ID_Limiter == 3){
            Interpolation2Face_Quick_K(i, j, k, n,
                                       f_Face,
                                       f,
                                       V,
                                       ID_Dimension);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            Interpolation2Face_WENO3_K(i, j, k, n,
                                       f_Face,
                                       f,
                                       V,
                                       ID_Dimension);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            Interpolation2Face_WENO3B_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        f_max,
                                        f_min,
                                        V,
                                        ID_Dimension,
                                        ID_BoundType);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            Interpolation2Face_WENO5_K(i, j, k, n,
                                       f_Face,
                                       f,
                                       V,
                                       ID_Dimension);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            Interpolation2Face_WENO5B_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        f_max,
                                        f_min,
                                        V,
                                        ID_Dimension,
                                        ID_Poly,
                                        ID_BoundType);
        }
    }
    
}

#endif