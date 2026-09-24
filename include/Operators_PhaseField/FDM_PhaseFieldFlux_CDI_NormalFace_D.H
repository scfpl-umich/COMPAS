// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_CDI_NormalFace_D_H_
#define FDM_PhaseFieldFlux_CDI_NormalFace_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_CDI_NormalFace_D(amrex::Real& PF_Flux,
                                         amrex::Real const& PF_c,
                                         amrex::Real const& PF_xm,
#if (AMREX_SPACEDIM > 1)
                                         amrex::Real const& PF_xp,
                                         amrex::Real const& PF_ym,
                                         amrex::Real const& PF_yp,
                                         amrex::Real const& PF_xm_ym,
                                         amrex::Real const& PF_xm_yp,
                                         amrex::Real const& PF_xp_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                                         amrex::Real const& PF_zm,
                                         amrex::Real const& PF_zp,
                                         amrex::Real const& PF_xm_zm,
                                         amrex::Real const& PF_xm_zp,
                                         amrex::Real const& PF_xp_zm,
                                         amrex::Real const& PF_ym_zm,
                                         amrex::Real const& PF_ym_zp,
                                         amrex::Real const& PF_yp_zm,
#endif
                                         amrex::Real const& Mob,
                                         amrex::Real const& Eta,
                            AMREX_D_DECL(amrex::Real const& dx,
                                         amrex::Real const& dy,
                                         amrex::Real const& dz),
                                         int ID_Dimension)
{
    amrex::Real Grad_PF;
    amrex::Real PF_m;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal_PF;
    //x
    if (ID_Dimension == 0){
        FDM_GradientX_XFaceFromCenter_D(Grad_PF,
                                        PF_c, PF_xm,
                                        dx);
        PF_m = PF_xm;
        
        FDM_Normal_XFaceFromCenter_D(Normal_PF,
                                     PF_c,
                                     PF_xm,
#if (AMREX_SPACEDIM > 1)
                                     PF_ym, PF_yp, PF_xm_ym, PF_xm_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                     PF_zm, PF_zp, PF_xm_zm, PF_xm_zp,
#endif
                        AMREX_D_DECL(dx, dy, dz));
    }
#if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        FDM_GradientY_YFaceFromCenter_D(Grad_PF,
                                        PF_c, PF_ym,
                                        dy);
        
        PF_m = PF_ym;
        
        FDM_Normal_YFaceFromCenter_D(Normal_PF,
                                     PF_c,
                                     PF_xm, PF_xp,
                                     PF_ym, PF_xm_ym, PF_xp_ym,
#if (AMREX_SPACEDIM > 2)
                                     PF_zm, PF_zp, PF_ym_zm, PF_ym_zp,
#endif
                        AMREX_D_DECL(dx, dy, dz));
    }
#endif
#if (AMREX_SPACEDIM > 2)
    //z
    else if (ID_Dimension == 2){
        FDM_GradientZ_ZFaceFromCenter_D(Grad_PF,
                                        PF_c, PF_zm,
                                        dz);
        
        PF_m = PF_zm;
        
        FDM_Normal_ZFaceFromCenter_D(Normal_PF,
                                     PF_c,
                                     PF_xm, PF_xp,
                                     PF_ym, PF_yp,
                                     PF_zm, PF_xm_zm, PF_xp_zm, PF_ym_zm, PF_yp_zm,
                        AMREX_D_DECL(dx, dy, dz));
    }
#endif
    
    amrex::Real PF;
    Interpolation2Face_Linear_D(PF,
                                PF_m, PF_c);
    
    amrex::Real Normal = Normal_PF[ID_Dimension];
    
    PhaseFieldFlux_CDI(PF_Flux,
                       Grad_PF, PF, Normal,
                       Mob, Eta);
}

#endif