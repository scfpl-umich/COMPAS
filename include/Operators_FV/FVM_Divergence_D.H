// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Divergence_D_H_
#define FVM_Divergence_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Divergence_D(amrex::Real& Div_f,
                      amrex::Real const& fx_fp,
                      amrex::Real const& fx_fm,
#if (AMREX_SPACEDIM > 1)
                      amrex::Real const& fy_fp,
                      amrex::Real const& fy_fm,
#endif
#if (AMREX_SPACEDIM > 2)
                      amrex::Real const& fz_fp,
                      amrex::Real const& fz_fm,
#endif
         AMREX_D_DECL(amrex::Real const& dx,
                      amrex::Real const& dy,
                      amrex::Real const& dz))
{
    Div_f = AMREX_D_TERM( (fx_fp - fx_fm)/dx,
                        + (fy_fp - fy_fm)/dy,
                        + (fz_fp - fz_fm)/dz );
}

#endif