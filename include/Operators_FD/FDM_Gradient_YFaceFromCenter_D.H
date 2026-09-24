// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#if (AMREX_SPACEDIM > 1)

#ifndef FDM_Gradient_YFaceFromCenter_D_H_
#define FDM_Gradient_YFaceFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_YFaceFromCenter_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Grad_f,
                                    amrex::Real const& f_c,
                                    amrex::Real const& f_xm,
                                    amrex::Real const& f_xp,
                                    amrex::Real const& f_ym,
                                    amrex::Real const& f_xm_ym,
                                    amrex::Real const& f_xp_ym,
#if (AMREX_SPACEDIM > 2)
                                    amrex::Real const& f_zm,
                                    amrex::Real const& f_zp,
                                    amrex::Real const& f_ym_zm,
                                    amrex::Real const& f_ym_zp,
#endif
                       AMREX_D_DECL(amrex::Real const& dx,
                                    amrex::Real const& dy,
                                    amrex::Real const& dz))
{
    amrex::Real Grad_f0;
    
    //dfdx(i,j-1/2,k,n)
    FDM_GradientX_YFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_xm,
                                    f_xp,
                                    f_ym,
                                    f_xm_ym,
                                    f_xp_ym,
                                    dx);
    Grad_f[0] = Grad_f0;
    
    //dfdy(i,j-1/2,k,n)
    FDM_GradientY_YFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_ym,
                                    dy);
    Grad_f[1] = Grad_f0;
    
#if (AMREX_SPACEDIM > 2)
    //dfdz(i,j-1/2,k,n)
    FDM_GradientZ_YFaceFromCenter_D(Grad_f0,
                                    f_c,
                                    f_ym,
                                    f_zm,
                                    f_zp,
                                    f_ym_zm,
                                    f_ym_zp,
                                    dz);
    Grad_f[2] = Grad_f0;
#endif
    
}

#endif

#endif