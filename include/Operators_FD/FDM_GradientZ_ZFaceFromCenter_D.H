// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_GradientZ_ZFaceFromCenter_D_H_
#define FDM_GradientZ_ZFaceFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_GradientZ_ZFaceFromCenter_D(amrex::Real& Grad_f,
                                     amrex::Real const& f_c,
                                     amrex::Real const& f_zm,
                                     amrex::Real const& dz)
{
    //dfdz(i,j,k-1/2,n)
    amrex::Real f_fm = f_zm; //f_fm = f(i,j,k-1,n)
    amrex::Real f_fp = f_c;  //f_fp = f(i,j,k,n)
    
    FDM_Gradient_D(Grad_f,
                   f_fm,
                   f_fp,
                   dz);
    
}

#endif