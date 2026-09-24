// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Normal_CenterFromCenter_K_H_
#define FDM_Normal_CenterFromCenter_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Normal_CenterFromCenter_K(int i, int j, int k, int n,
                                   amrex::Array<amrex::Real,AMREX_SPACEDIM>& Normal_f,
                                   amrex::Array4<amrex::Real const> const& f,
                      AMREX_D_DECL(amrex::Real const& dx,
                                   amrex::Real const& dy,
                                   amrex::Real const& dz))
{
    amrex::Real f_c  = f(i,j,k,n);
    amrex::Real f_xm = f(i-1,j,k,n);
    amrex::Real f_xp = f(i+1,j,k,n);
#if (AMREX_SPACEDIM > 1)
    amrex::Real f_ym = f(i,j-1,k,n);
    amrex::Real f_yp = f(i,j+1,k,n);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real f_zm = f(i,j,k-1,n);
    amrex::Real f_zp = f(i,j,k+1,n);
#endif
    FDM_Normal_CenterFromCenter_D(Normal_f,
                                  f_c,
                                  f_xm, f_xp,
#if (AMREX_SPACEDIM > 1)
                                  f_ym, f_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                  f_zm, f_zp,
#endif
                     AMREX_D_DECL(dx, dy, dz));
    
}

#endif