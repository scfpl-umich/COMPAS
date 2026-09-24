// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Divergence_K_H_
#define FVM_Divergence_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Divergence_K(int i, int j, int k, int n,
                      amrex::Real& Div_f,
                      AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fx,
                                   amrex::Array4<amrex::Real const> const& fy,
                                   amrex::Array4<amrex::Real const> const& fz),
                      AMREX_D_DECL(amrex::Real const& dx,
                                   amrex::Real const& dy,
                                   amrex::Real const& dz))
{
    amrex::Real fx_fp = fx(i+1,j,k,n);
    amrex::Real fx_fm = fx(i,j,k,n);
#if (AMREX_SPACEDIM > 1)
    amrex::Real fy_fp = fy(i,j+1,k,n);
    amrex::Real fy_fm = fy(i,j,k,n);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real fz_fp = fz(i,j,k+1,n);
    amrex::Real fz_fm = fz(i,j,k,n);
#endif
    
    FVM_Divergence_D(Div_f,
                     fx_fp, fx_fm,
#if (AMREX_SPACEDIM > 1)
                     fy_fp, fy_fm,
#endif
#if (AMREX_SPACEDIM > 2)
                     fz_fp, fz_fm,
#endif
         AMREX_D_DECL(dx, dy, dz));
    
}

#endif