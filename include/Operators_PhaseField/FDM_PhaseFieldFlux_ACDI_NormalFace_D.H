// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_ACDI_NormalFace_D_H_
#define FDM_PhaseFieldFlux_ACDI_NormalFace_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_ACDI_NormalFace_D(amrex::Real& PF_Flux,
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
    int ID_Derivative = 0;
    amrex::Real Psi_c;
    Mapping_Backward(Psi_c,
                     PF_c, Eta, ID_Derivative);
    amrex::Real Psi_xm;
    Mapping_Backward(Psi_xm,
                     PF_xm, Eta, ID_Derivative);
#if (AMREX_SPACEDIM > 1)
    amrex::Real Psi_xp;
    Mapping_Backward(Psi_xp,
                     PF_xp, Eta, ID_Derivative);
    amrex::Real Psi_ym;
    Mapping_Backward(Psi_ym,
                     PF_ym, Eta, ID_Derivative);
    amrex::Real Psi_yp;
    Mapping_Backward(Psi_yp,
                     PF_yp, Eta, ID_Derivative);
    amrex::Real Psi_xm_ym;
    Mapping_Backward(Psi_xm_ym,
                     PF_xm_ym, Eta, ID_Derivative);
    amrex::Real Psi_xm_yp;
    Mapping_Backward(Psi_xm_yp,
                     PF_xm_yp, Eta, ID_Derivative);
    amrex::Real Psi_xp_ym;
    Mapping_Backward(Psi_xp_ym,
                     PF_xp_ym, Eta, ID_Derivative);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Psi_zm;
    Mapping_Backward(Psi_zm,
                     PF_zm, Eta, ID_Derivative);
    amrex::Real Psi_zp;
    Mapping_Backward(Psi_zp,
                     PF_zp, Eta, ID_Derivative);
    amrex::Real Psi_xm_zm;
    Mapping_Backward(Psi_xm_zm,
                     PF_xm_zm, Eta, ID_Derivative);
    amrex::Real Psi_xm_zp;
    Mapping_Backward(Psi_xm_zp,
                     PF_xm_zp, Eta, ID_Derivative);
    amrex::Real Psi_xp_zm;
    Mapping_Backward(Psi_xp_zm,
                     PF_xp_zm, Eta, ID_Derivative);
    amrex::Real Psi_ym_zm;
    Mapping_Backward(Psi_ym_zm,
                     PF_ym_zm, Eta, ID_Derivative);
    amrex::Real Psi_ym_zp;
    Mapping_Backward(Psi_ym_zp,
                     PF_ym_zp, Eta, ID_Derivative);
    amrex::Real Psi_yp_zm;
    Mapping_Backward(Psi_yp_zm,
                     PF_yp_zm, Eta, ID_Derivative);
#endif
    
    amrex::Real Grad_PF;
    amrex::Real Psi_m;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal_Psi;
    //x
    if (ID_Dimension == 0){
        FDM_GradientX_XFaceFromCenter_D(Grad_PF,
                                        PF_c, PF_xm,
                                        dx);
        Psi_m = Psi_xm;
        
        FDM_Normal_XFaceFromCenter_D(Normal_Psi,
                                     Psi_c,
                                     Psi_xm,
#if (AMREX_SPACEDIM > 1)
                                     Psi_ym, Psi_yp, Psi_xm_ym, Psi_xm_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                     Psi_zm, Psi_zp, Psi_xm_zm, Psi_xm_zp,
#endif
                        AMREX_D_DECL(dx, dy, dz));
    }
#if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        FDM_GradientY_YFaceFromCenter_D(Grad_PF,
                                        PF_c, PF_ym,
                                        dy);
        
        Psi_m = Psi_ym;
        
        FDM_Normal_YFaceFromCenter_D(Normal_Psi,
                                     Psi_c,
                                     Psi_xm, Psi_xp,
                                     Psi_ym, Psi_xm_ym, Psi_xp_ym,
#if (AMREX_SPACEDIM > 2)
                                     Psi_zm, Psi_zp, Psi_ym_zm, Psi_ym_zp,
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
        
        Psi_m = Psi_zm;
        
        FDM_Normal_ZFaceFromCenter_D(Normal_Psi,
                                     Psi_c,
                                     Psi_xm, Psi_xp,
                                     Psi_ym, Psi_yp,
                                     Psi_zm, Psi_xm_zm, Psi_xp_zm, Psi_ym_zm, Psi_yp_zm,
                        AMREX_D_DECL(dx, dy, dz));
    }
#endif
    amrex::Real Psi;
    Interpolation2Face_Linear_D(Psi,
                                Psi_m, Psi_c);
    
    amrex::Real Normal = Normal_Psi[ID_Dimension];
    
    PhaseFieldFlux_ACDI(PF_Flux,
                        Grad_PF, Psi, Normal,
                        Mob, Eta);
}

#endif