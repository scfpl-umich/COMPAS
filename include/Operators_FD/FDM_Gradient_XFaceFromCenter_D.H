// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Gradient_XFaceFromCenter_D_H_
#define FDM_Gradient_XFaceFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_XFaceFromCenter_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Grad_f,
                                    amrex::Real const& f_c,
                                    amrex::Real const& f_xm,
#if (AMREX_SPACEDIM > 1)
                                    amrex::Real const& f_ym,
                                    amrex::Real const& f_yp,
                                    amrex::Real const& f_xm_ym,
                                    amrex::Real const& f_xm_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                    amrex::Real const& f_zm,
                                    amrex::Real const& f_zp,
                                    amrex::Real const& f_xm_zm,
                                    amrex::Real const& f_xm_zp,
#endif
                       AMREX_D_DECL(amrex::Real const& dx,
                                    amrex::Real const& dy,
                                    amrex::Real const& dz))
{
    amrex::Real Grad_f0;
    
    //dfdx(i-1/2,j,k,n)
    FDM_GradientX_XFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_xm,
                                    dx);
    Grad_f[0] = Grad_f0;
    
#if (AMREX_SPACEDIM > 1)
    //dfdy(i-1/2,j,k,n)
    FDM_GradientY_XFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_xm,
                                    f_ym,
                                    f_yp,
                                    f_xm_ym,
                                    f_xm_yp,
                                    dy);
    Grad_f[1] = Grad_f0;
#endif
    
#if (AMREX_SPACEDIM > 2)
    //dfdz(i-1/2,j,k,n)
    FDM_GradientZ_XFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_xm,
                                    f_zm,
                                    f_zp,
                                    f_xm_zm,
                                    f_xm_zp,
                                    dz);
    Grad_f[2] = Grad_f0;
#endif
    
}

#endif