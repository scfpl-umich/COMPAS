// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Gradient_FaceFromCenter_K_H_
#define FDM_Gradient_FaceFromCenter_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_FaceFromCenter_K(int i, int j, int k, int n,
                                   amrex::Real& Grad_f,
                                   amrex::Array4<amrex::Real const> const& f,
                                   amrex::Real const& dx,
                                   int ID_Dimension)
{
    amrex::Real f_m;
    amrex::Real f_c;
    
    f_c = f(i,j,k,n);
    //x
    if (ID_Dimension == 0){
        //dfdx(i-1/2,j,k,n)
        f_m = f(i-1,j,k,n); //f_xm
        FDM_GradientX_XFaceFromCenter_D(Grad_f,
                                        f_c,
                                        f_m,
                                        dx);
    }
    //y
    else if (ID_Dimension == 1){
        //dfdy(i,j-1/2,k,n)
        f_m = f(i,j-1,k,n); //f_ym
        FDM_GradientY_YFaceFromCenter_D(Grad_f,
                                        f_c,
                                        f_m,
                                        dx);
    }
    //z
    else if (ID_Dimension == 2){
        //dfdz(i,j,k-1/2,n)
        f_m = f(i,j,k-1,n); //f_zm
        FDM_GradientZ_ZFaceFromCenter_D(Grad_f,
                                        f_c,
                                        f_m,
                                        dx);
    }
    
}

#endif