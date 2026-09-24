// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Normal_CenterFromCenter_D_H_
#define FDM_Normal_CenterFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Normal_CenterFromCenter_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Normal_f,
                                   amrex::Real const& f_c,
                                   amrex::Real const& f_xm,
                                   amrex::Real const& f_xp,
#if (AMREX_SPACEDIM > 1)
                                   amrex::Real const& f_ym,
                                   amrex::Real const& f_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                   amrex::Real const& f_zm,
                                   amrex::Real const& f_zp,
#endif
                      AMREX_D_DECL(amrex::Real const& dx,
                                   amrex::Real const& dy,
                                   amrex::Real const& dz))
{
    amrex::Real Sum = 0.0;
    amrex::Real Grad_f;
    //dfdx
    FDM_Gradient_CenterFromCenter_D(Grad_f,
                                    f_xm, f_c, f_xp,
                                    dx);
    Sum += Grad_f*Grad_f;
    Normal_f[0] = Grad_f;
#if (AMREX_SPACEDIM > 1)
    //dfdy
    FDM_Gradient_CenterFromCenter_D(Grad_f,
                                    f_ym, f_c, f_yp,
                                    dy);
    Sum += Grad_f*Grad_f;
    Normal_f[1] = Grad_f;
#endif
#if (AMREX_SPACEDIM > 2)
    //dfdz
    FDM_Gradient_CenterFromCenter_D(Grad_f,
                                    f_zm, f_c, f_zp,
                                    dz);
    Sum += Grad_f*Grad_f;
    Normal_f[2] = Grad_f;
#endif
    if (Sum == 0.0){
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
            Normal_f[idim] = 0.0;
        }
    }
    else{
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
            Normal_f[idim] /= std::sqrt(Sum);
        }
    }
    
}

#endif