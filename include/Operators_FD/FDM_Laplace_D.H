// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Laplace_D_H_
#define FDM_Laplace_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Laplace_D(amrex::Real& Laplace_f,
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
    amrex::Real Gradx_fp;
    FDM_GradientX_XFaceFromCenter_D(Gradx_fp,
                                    f_xp, f_c,
                                    dx);
    amrex::Real Gradx_fm;
    FDM_GradientX_XFaceFromCenter_D(Gradx_fm,
                                    f_c, f_xm,
                                    dx);
#if (AMREX_SPACEDIM > 1)
    amrex::Real Grady_fp;
    FDM_GradientY_YFaceFromCenter_D(Grady_fp,
                                    f_yp, f_c,
                                    dy);
    amrex::Real Grady_fm;
    FDM_GradientY_YFaceFromCenter_D(Grady_fm,
                                    f_c, f_ym,
                                    dy);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Gradz_fp;
    FDM_GradientZ_ZFaceFromCenter_D(Gradz_fp,
                                    f_zp, f_c,
                                    dz);
    amrex::Real Gradz_fm;
    FDM_GradientZ_ZFaceFromCenter_D(Gradz_fm,
                                    f_c, f_zm,
                                    dz);
#endif
    FVM_Divergence_D(Laplace_f,
                     Gradx_fp, Gradx_fm,
#if (AMREX_SPACEDIM > 1)
                     Grady_fp, Grady_fm,
#endif
#if (AMREX_SPACEDIM > 2)
                     Gradz_fp, Gradz_fm,
#endif
        AMREX_D_DECL(dx, dy, dz));
    
}

#endif