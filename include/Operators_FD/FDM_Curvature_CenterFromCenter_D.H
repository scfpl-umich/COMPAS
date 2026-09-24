// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Curvature_CenterFromCenter_D_H_
#define FDM_Curvature_CenterFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Curvature_CenterFromCenter_D(amrex::Real& Curvature_f,
                                      amrex::Real const& f_c,
                                      amrex::Real const& f_xm,
                                      amrex::Real const& f_xp,
    #if (AMREX_SPACEDIM > 1)
                                      amrex::Real const& f_ym,
                                      amrex::Real const& f_yp,
                                      amrex::Real const& f_xm_ym,
                                      amrex::Real const& f_xm_yp,
                                      amrex::Real const& f_xp_ym,
                                      amrex::Real const& f_xp_yp,
    #endif
    #if (AMREX_SPACEDIM > 2)
                                      amrex::Real const& f_zm,
                                      amrex::Real const& f_zp,
                                      amrex::Real const& f_xm_zm,
                                      amrex::Real const& f_xm_zp,
                                      amrex::Real const& f_xp_zm,
                                      amrex::Real const& f_xp_zp,
                                      amrex::Real const& f_ym_zm,
                                      amrex::Real const& f_ym_zp,
                                      amrex::Real const& f_yp_zm,
                                      amrex::Real const& f_yp_zp,
    #endif
                         AMREX_D_DECL(amrex::Real const& dx,
                                      amrex::Real const& dy,
                                      amrex::Real const& dz))
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal_f;

    //nx_fm
    FDM_Normal_XFaceFromCenter_D(Normal_f,
                                 f_c,
                                 f_xm,
#if (AMREX_SPACEDIM > 1)
                                 f_ym, f_yp, f_xm_ym, f_xm_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                 f_zm, f_zp, f_xm_zm, f_xm_zp,
#endif
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real nx_fm = Normal_f[0];

    //nx_fp
    FDM_Normal_XFaceFromCenter_D(Normal_f,
                                 f_xp,
                                 f_c,
#if (AMREX_SPACEDIM > 1)
                                 f_xp_ym, f_xp_yp, f_ym, f_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                 f_xp_zm, f_xp_zp, f_zm, f_zp,
#endif
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real nx_fp = Normal_f[0];

#if (AMREX_SPACEDIM > 1)
    //ny_fm
    FDM_Normal_YFaceFromCenter_D(Normal_f,
                                 f_c,
                                 f_xm, f_xp,
                                 f_ym, f_xm_ym, f_xp_ym,
#if (AMREX_SPACEDIM > 2)
                                 f_zm, f_zp, f_ym_zm, f_ym_zp,
#endif
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real ny_fm = Normal_f[1];

    //ny_fp
    FDM_Normal_YFaceFromCenter_D(Normal_f,
                                 f_yp,
                                 f_xm_yp, f_xp_yp,
                                 f_c, f_xm, f_xp,
#if (AMREX_SPACEDIM > 2)
                                 f_yp_zm, f_yp_zp, f_zm, f_zp,
#endif
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real ny_fp = Normal_f[1];
#endif

#if (AMREX_SPACEDIM > 2)
    //nz_fm
    FDM_Normal_ZFaceFromCenter_D(Normal_f,
                                 f_c,
                                 f_xm, f_xp,
                                 f_ym, f_yp,
                                 f_zm, f_xm_zm, f_xp_zm, f_ym_zm, f_yp_zm,
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real nz_fm = Normal_f[2];
    
    //nz_fp
    FDM_Normal_ZFaceFromCenter_D(Normal_f,
                                 f_zp,
                                 f_xm_zp, f_xp_zp,
                                 f_ym_zp, f_yp_zp,
                                 f_c, f_xm, f_xp, f_ym, f_yp,
                    AMREX_D_DECL(dx, dy, dz));
    amrex::Real nz_fp = Normal_f[2];
#endif

    amrex::Real Div_Normal;
    FVM_Divergence_D(Div_Normal,
                     nx_fp, nx_fm,
#if (AMREX_SPACEDIM > 1)
                     ny_fp, ny_fm,
#endif
#if (AMREX_SPACEDIM > 2)
                     nz_fp, nz_fm,
#endif
        AMREX_D_DECL(dx, dy, dz));

    Curvature_f = -Div_Normal;
    
}

#endif