// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_D_H_
#define Reconstruction_D_H_

template<typename T_Stencil>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_D(amrex::Real& f_Reconstruction,
                      amrex::Real const& x,
                      T_Stencil const& f_Stencil,
                      amrex::Real const& f_max,
                      amrex::Real const& f_min,
                      amrex::Real const& x_c,
                      amrex::Real const& dx,
                      int ID_WENO,
                      int ID_Limiter,
                      int ID_Poly,
                      int ID_BoundType)
{
    //non-WENO
    if (ID_WENO == 0){
        //Godunov (1st order)
        if (ID_Limiter == 0){
            amrex::Real f_c = f_Stencil[0];
            Reconstruction_Godunov_D(f_Reconstruction,
                                     x,
                                     f_c);
        }
        //MUSCL (2nd order)
        else{
            amrex::Real f_m = f_Stencil[0];
            amrex::Real f_c = f_Stencil[1];
            amrex::Real f_p = f_Stencil[2];
            Reconstruction_MUSCL_D(f_Reconstruction,
                                   x,
                                   f_m, f_c, f_p,
                                   x_c, dx,
                                   ID_Limiter);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            amrex::Real f_m = f_Stencil[0];
            amrex::Real f_c = f_Stencil[1];
            amrex::Real f_p = f_Stencil[2];
            Reconstruction_WENO3_D(f_Reconstruction,
                                   x,
                                   f_m, f_c, f_p,
                                   x_c, dx);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            amrex::Real f_m = f_Stencil[0];
            amrex::Real f_c = f_Stencil[1];
            amrex::Real f_p = f_Stencil[2];
            Reconstruction_WENO3B_D(f_Reconstruction,
                                    x,
                                    f_m, f_c, f_p,
                                    f_max, f_min,
                                    x_c, dx,
                                    ID_BoundType);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            amrex::Real f_mm = f_Stencil[0];
            amrex::Real f_m  = f_Stencil[1];
            amrex::Real f_c  = f_Stencil[2];
            amrex::Real f_p  = f_Stencil[3];
            amrex::Real f_pp = f_Stencil[4];
            Reconstruction_WENO5_D(f_Reconstruction,
                                   x,
                                   f_mm, f_m, f_c, f_p, f_pp,
                                   x_c, dx,
                                   ID_Poly);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            amrex::Real f_mm = f_Stencil[0];
            amrex::Real f_m  = f_Stencil[1];
            amrex::Real f_c  = f_Stencil[2];
            amrex::Real f_p  = f_Stencil[3];
            amrex::Real f_pp = f_Stencil[4];
            Reconstruction_WENO5B_D(f_Reconstruction,
                                    x,
                                    f_mm, f_m, f_c, f_p, f_pp,
                                    f_max, f_min,
                                    x_c, dx,
                                    ID_Poly, ID_BoundType);
        }
        
    }
    
}

#endif