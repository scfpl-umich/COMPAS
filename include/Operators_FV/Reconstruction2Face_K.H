// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_K_H_
#define Reconstruction2Face_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_K(int i, int j, int k, int n,
                           amrex::Real& f_Face,
                           amrex::Array4<amrex::Real const> const& f,
                           amrex::Real const& f_max,
                           amrex::Real const& f_min,
                           amrex::Real const& dx,
                           int ID_Dimension,
                           int ID_Face,
                           int ID_WENO,
                           int ID_Limiter,
                           int ID_Poly,
                           int ID_BoundType)
{
    //non-WENO
    if (ID_WENO == 0){
        //Godunov (1st order)
        if (ID_Limiter == 0){
            Reconstruction2Face_Godunov_K(i, j, k, n,
                                          f_Face,
                                          f);
        }
        //MUSCL (2nd order)
        else{
            Reconstruction2Face_MUSCL_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        dx,
                                        ID_Dimension,
                                        ID_Face,
                                        ID_Limiter);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            Reconstruction2Face_WENO3_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        ID_Dimension,
                                        ID_Face);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            Reconstruction2Face_WENO3B_K(i, j, k, n,
                                         f_Face,
                                         f,
                                         f_max,
                                         f_min,
                                         ID_Dimension,
                                         ID_Face,
                                         ID_BoundType);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            Reconstruction2Face_WENO5_K(i, j, k, n,
                                        f_Face,
                                        f,
                                        ID_Dimension,
                                        ID_Face);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            Reconstruction2Face_WENO5B_K(i, j, k, n,
                                         f_Face,
                                         f,
                                         f_max,
                                         f_min,
                                         ID_Dimension,
                                         ID_Face,
                                         ID_Poly,
                                         ID_BoundType);
        }
    }
    
}

#endif