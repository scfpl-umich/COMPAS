// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#if (AMREX_SPACEDIM > 2)

#ifndef FDM_Gradient_ZFaceFromCenter_D_H_
#define FDM_Gradient_ZFaceFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_ZFaceFromCenter_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Grad_f,
                                    amrex::Real const& f_c,
                                    amrex::Real const& f_xm,
                                    amrex::Real const& f_xp,
                                    amrex::Real const& f_ym,
                                    amrex::Real const& f_yp,
                                    amrex::Real const& f_zm,
                                    amrex::Real const& f_xm_zm,
                                    amrex::Real const& f_xp_zm,
                                    amrex::Real const& f_ym_zm,
                                    amrex::Real const& f_yp_zm,
                       AMREX_D_DECL(amrex::Real const& dx,
                                    amrex::Real const& dy,
                                    amrex::Real const& dz))
{
    amrex::Real Grad_f0;
    
    //dfdx(i,j,k-1/2,n)
    FDM_GradientX_ZFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_xm,
                                    f_xp,
                                    f_zm,
                                    f_xm_zm,
                                    f_xp_zm,
                                    dx);
    Grad_f[0] = Grad_f0;
    
    //dfdy(i,j,k-1/2,n)
    FDM_GradientY_ZFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_ym,
                                    f_yp,
                                    f_zm,
                                    f_ym_zm,
                                    f_yp_zm,
                                    dy);
    Grad_f[1] = Grad_f0;
    
    //dfdz(i,j,k-1/2,n)
    FDM_GradientZ_ZFaceFromCenter_D(Grad_f0,
                                   f_c,
                                   f_zm,
                                   dz);
    Grad_f[2] = Grad_f0;
    
}

#endif

#endif