// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_K_H_
#define Reconstruction_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_K(int i, int j, int k, int n,
                      amrex::Real& f_Reconstruction,
                      amrex::Real const& x,
                      amrex::Array4<amrex::Real const> const& f,
                      amrex::Real const& f_max,
                      amrex::Real const& f_min,
                      amrex::Real const& dx,
                      amrex::Real const& prob_lo,
                      int ID_Dimension,
                      int ID_WENO,
                      int ID_Limiter,
                      int ID_Poly,
                      int ID_BoundType)
{
    //non-WENO
    if (ID_WENO == 0){
        //Godunov (1st order)
        if (ID_Limiter == 0){
            Reconstruction_Godunov_K(i, j, k, n,
                                     f_Reconstruction,
                                     x,
                                     f);
        }
        //MUSCL (2nd order)
        else{
            Reconstruction_MUSCL_K(i, j, k, n,
                                   f_Reconstruction,
                                   x,
                                   f,
                                   dx,
                                   prob_lo,
                                   ID_Dimension,
                                   ID_Limiter);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            Reconstruction_WENO3_K(i, j, k, n,
                                   f_Reconstruction,
                                   x,
                                   f,
                                   dx,
                                   prob_lo,
                                   ID_Dimension);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            Reconstruction_WENO3B_K(i, j, k, n,
                                    f_Reconstruction,
                                    x,
                                    f,
                                    f_max,
                                    f_min,
                                    dx,
                                    prob_lo,
                                    ID_Dimension,
                                    ID_BoundType);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            Reconstruction_WENO5_K(i, j, k, n,
                                   f_Reconstruction,
                                   x,
                                   f,
                                   dx,
                                   prob_lo,
                                   ID_Dimension,
                                   ID_Poly);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            Reconstruction_WENO5B_K(i, j, k, n,
                                    f_Reconstruction,
                                    x,
                                    f,
                                    f_max,
                                    f_min,
                                    dx,
                                    prob_lo,
                                    ID_Dimension,
                                    ID_Poly,
                                    ID_BoundType);
        }
        
    }
    
}

#endif