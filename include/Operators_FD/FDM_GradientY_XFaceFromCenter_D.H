// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_GradientY_XFaceFromCenter_D_H_
#define FDM_GradientY_XFaceFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_GradientY_XFaceFromCenter_D(amrex::Real& Grad_f,
                                     amrex::Real const& f_c,
                                     amrex::Real const& f_xm,
                                     amrex::Real const& f_ym,
                                     amrex::Real const& f_yp,
                                     amrex::Real const& f_xm_ym,
                                     amrex::Real const& f_xm_yp,
                                     amrex::Real const& dy)
{
    //dfdy(i-1/2,j,k,n)
    amrex::Real f_fm; //f_fm = f(i-1/2,j-1/2,k,n)
    Interpolation2Edge_Linear_D(f_fm,
                                f_c,
                                f_xm,
                                f_ym,
                                f_xm_ym);
    
    amrex::Real f_fp; //f_fp = f(i-1/2,j+1/2,k,n)
    Interpolation2Edge_Linear_D(f_fp,
                                f_c,
                                f_xm,
                                f_yp,
                                f_xm_yp);
    
    FDM_Gradient_D(Grad_f,
                   f_fm,
                   f_fp,
                   dy);
    
}

#endif