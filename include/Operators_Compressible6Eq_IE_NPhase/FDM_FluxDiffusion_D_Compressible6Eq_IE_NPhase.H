// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_
#define FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE

void FDM_FluxDiffusion_D(amrex::Array<amrex::Real, NSTATE>& Flux,
                         amrex::Real& dmax,
                         amrex::Array<amrex::Real,AMREX_SPACEDIM>& tau,
                         amrex::Real const& u_c,
                         amrex::Real const& u_xm,
                         amrex::Real const& T_c,
                         amrex::Real const& T_xm,
                         amrex::Real const& mu_c,
                         amrex::Real const& mu_xm,
                         amrex::Real const& mub_c,
                         amrex::Real const& mub_xm,
                         amrex::Real const& Kappa_c,
                         amrex::Real const& Kappa_xm,
#if (AMREX_SPACEDIM > 1)
                         amrex::Real const& u_xp, amrex::Real const& u_ym, amrex::Real const& u_yp,
                         amrex::Real const& u_xm_ym, amrex::Real const& u_xm_yp, amrex::Real const& u_xp_ym,
                         amrex::Real const& v_c,
                         amrex::Real const& v_xm, amrex::Real const& v_xp, amrex::Real const& v_ym, amrex::Real const& v_yp,
                         amrex::Real const& v_xm_ym, amrex::Real const& v_xm_yp, amrex::Real const& v_xp_ym,
                         amrex::Real const& T_ym,
                         amrex::Real const& mu_ym,
                         amrex::Real const& mub_ym,
                         amrex::Real const& Kappa_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                         amrex::Real const& u_zm, amrex::Real const& u_zp,
                         amrex::Real const& u_xm_zm, amrex::Real const& u_xm_zp, amrex::Real const& u_xp_zm,
                         amrex::Real const& v_zm, amrex::Real const& v_zp,
                         amrex::Real const& v_ym_zm, amrex::Real const& v_ym_zp, amrex::Real const& v_yp_zm,
                         amrex::Real const& w_c,
                         amrex::Real const& w_xm, amrex::Real const& w_xp, amrex::Real const& w_ym, amrex::Real const& w_yp, amrex::Real const& w_zm, amrex::Real const& w_zp,
                         amrex::Real const& w_xm_zm, amrex::Real const& w_xm_zp, amrex::Real const& w_xp_zm,
                         amrex::Real const& w_ym_zm, amrex::Real const& w_ym_zp, amrex::Real const& w_yp_zm,
                         amrex::Real const& T_zm,
                         amrex::Real const& mu_zm,
                         amrex::Real const& mub_zm,
                         amrex::Real const& Kappa_zm,
#endif
            AMREX_D_DECL(amrex::Real const& dx,
                         amrex::Real const& dy,
                         amrex::Real const& dz),
                         int ID_Dimension)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Grad_Vel;
    amrex::Real Grad_T;
    amrex::Real mu;
    amrex::Real mub;
    amrex::Real Kappa;
/*
    u = Vel[0];
    du_dx = Grad_Vel[0][0];
#if (AMREX_SPACEDIM > 1)
    v = Vel[1];
    du_dy = Grad_Vel[1][0];
    dv_dx = Grad_Vel[0][1];
    dv_dy = Grad_Vel[1][1];
#endif
#if (AMREX_SPACEDIM > 2)
    w = Vel[2];
    du_dz = Grad_Vel[2][0];
    dv_dz = Grad_Vel[2][1];
    dw_dx = Grad_Vel[0][2];
    dw_dy = Grad_Vel[1][2];
    dw_dz = Grad_Vel[2][2];
#endif
*/
    //x
    if (ID_Dimension == 0){
        //Velocity and its gradient
        amrex::Real u;
        Interpolation2Face_Linear_D(u,
                                    u_xm,
                                    u_c);
        Vel[0] = u;
        
        amrex::Real du_dx;
        FDM_GradientX_XFaceFromCenter_D(du_dx,
                                        u_c,
                                        u_xm,
                                        dx);
        Grad_Vel[0][0] = du_dx;
#if (AMREX_SPACEDIM > 1)
        amrex::Real v;
        Interpolation2Face_Linear_D(v,
                                    v_xm,
                                    v_c);
        Vel[1] = v;
        
        amrex::Real du_dy;
        FDM_GradientY_XFaceFromCenter_D(du_dy,
                                        u_c,
                                        u_xm, u_ym, u_yp,
                                        u_xm_ym, u_xm_yp,
                                        dy);
        Grad_Vel[1][0] = du_dy;
        
        amrex::Real dv_dx;
        FDM_GradientX_XFaceFromCenter_D(dv_dx,
                                        v_c,
                                        v_xm,
                                        dx);
        Grad_Vel[0][1] = dv_dx;
        
        amrex::Real dv_dy;
        FDM_GradientY_XFaceFromCenter_D(dv_dy,
                                        v_c,
                                        v_xm, v_ym, v_yp,
                                        v_xm_ym, v_xm_yp,
                                        dy);
        Grad_Vel[1][1] = dv_dy;
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real w;
        Interpolation2Face_Linear_D(w,
                                    w_xm,
                                    w_c);
        Vel[2] = w;
        
        amrex::Real du_dz;
        FDM_GradientZ_XFaceFromCenter_D(du_dz,
                                        u_c,
                                        u_xm, u_zm, u_zp,
                                        u_xm_zm, u_xm_zp,
                                        dz);
        Grad_Vel[2][0] = du_dz;
        
        amrex::Real dw_dx;
        FDM_GradientX_XFaceFromCenter_D(dw_dx,
                                        w_c,
                                        w_xm,
                                        dx);
        Grad_Vel[0][2] = dw_dx;
        
        amrex::Real dw_dz;
        FDM_GradientZ_XFaceFromCenter_D(dw_dz,
                                        w_c,
                                        w_xm, w_zm, w_zp,
                                        w_xm_zm, w_xm_zp,
                                        dz);
        Grad_Vel[2][2] = dw_dz;
#endif
        //Temperature gradient
        amrex::Real dT_dx;
        FDM_GradientX_XFaceFromCenter_D(dT_dx,
                                        T_c,
                                        T_xm,
                                        dx);
        Grad_T = dT_dx;
        //Material property
        //mu
        Interpolation2Face_Linear_D(mu,
                                    mu_xm,
                                    mu_c);
        //mub
        Interpolation2Face_Linear_D(mub,
                                    mub_xm,
                                    mub_c);
        //Kappa
        Interpolation2Face_Linear_D(Kappa,
                                    Kappa_xm,
                                    Kappa_c);
    }
#if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        //Velocity and its gradient
        amrex::Real u;
        Interpolation2Face_Linear_D(u,
                                    u_ym,
                                    u_c);
        Vel[0] = u;
        
        amrex::Real v;
        Interpolation2Face_Linear_D(v,
                                    v_ym,
                                    v_c);
        Vel[1] = v;
        
        amrex::Real dv_dx;
        FDM_GradientX_YFaceFromCenter_D(dv_dx,
                                        v_c,
                                        v_xm, v_xp, v_ym,
                                        v_xm_ym, v_xp_ym,
                                        dx);
        Grad_Vel[0][1] = dv_dx;
        
        amrex::Real du_dx;
        FDM_GradientX_YFaceFromCenter_D(du_dx,
                                        u_c,
                                        u_xm, u_xp, u_ym,
                                        u_xm_ym, u_xp_ym,
                                        dx);
        Grad_Vel[0][0] = du_dx;
        
        amrex::Real du_dy;
        FDM_GradientY_YFaceFromCenter_D(du_dy,
                                        u_c,
                                        u_ym,
                                        dy);
        Grad_Vel[1][0] = du_dy;
        
        amrex::Real dv_dy;
        FDM_GradientY_YFaceFromCenter_D(dv_dy,
                                        v_c,
                                        v_ym,
                                        dy);
        Grad_Vel[1][1] = dv_dy;
#if (AMREX_SPACEDIM > 2)
        amrex::Real w;
        Interpolation2Face_Linear_D(w,
                                    w_ym,
                                    w_c);
        Vel[2] = w;
        
        amrex::Real dv_dz;
        FDM_GradientZ_YFaceFromCenter_D(dv_dz,
                                        v_c,
                                        v_ym, v_zm, v_zp,
                                        v_ym_zm, v_ym_zp,
                                        dz);
        Grad_Vel[2][1] = dv_dz;
        
        amrex::Real dw_dy;
        FDM_GradientY_YFaceFromCenter_D(dw_dy,
                                        w_c,
                                        w_ym,
                                        dy);
        Grad_Vel[1][2] = dw_dy;
        
        amrex::Real dw_dz;
        FDM_GradientZ_YFaceFromCenter_D(dw_dz,
                                        w_c,
                                        w_ym, w_zm, w_zp,
                                        w_ym_zm, w_ym_zp,
                                        dz);
        Grad_Vel[2][2] = dw_dz;
#endif
        //Temperature gradient
        amrex::Real dT_dy;
        FDM_GradientY_YFaceFromCenter_D(dT_dy,
                                        T_c,
                                        T_ym,
                                        dy);
        Grad_T = dT_dy;
        //Material property
        //mu
        Interpolation2Face_Linear_D(mu,
                                    mu_ym,
                                    mu_c);
        //mub
        Interpolation2Face_Linear_D(mub,
                                    mub_ym,
                                    mub_c);
        //Kappa
        Interpolation2Face_Linear_D(Kappa,
                                    Kappa_ym,
                                    Kappa_c);
    }
#endif
#if (AMREX_SPACEDIM > 2)
    //z
    else if (ID_Dimension == 2){
        //Velocity and its gradient
        amrex::Real u;
        Interpolation2Face_Linear_D(u,
                                    u_zm,
                                    u_c);
        Vel[0] = u;
        
        amrex::Real v;
        Interpolation2Face_Linear_D(v,
                                    v_zm,
                                    v_c);
        Vel[1] = v;
        
        amrex::Real w;
        Interpolation2Face_Linear_D(w,
                                    w_zm,
                                    w_c);
        Vel[2] = w;
        
        amrex::Real dw_dx;
        FDM_GradientX_ZFaceFromCenter_D(dw_dx,
                                        w_c,
                                        w_xm, w_xp, w_zm,
                                        w_xm_zm, w_xp_zm,
                                        dx);
        Grad_Vel[0][2] = dw_dx;
        
        amrex::Real du_dx;
        FDM_GradientX_ZFaceFromCenter_D(du_dx,
                                        u_c,
                                        u_xm, u_xp, u_zm,
                                        u_xm_zm, u_xp_zm,
                                        dx);
        Grad_Vel[0][0] = du_dx;
        
        amrex::Real du_dz;
        FDM_GradientZ_ZFaceFromCenter_D(du_dz,
                                        u_c,
                                        u_zm,
                                        dz);
        Grad_Vel[2][0] = du_dz;
        
        amrex::Real dw_dy;
        FDM_GradientY_ZFaceFromCenter_D(dw_dy,
                                        w_c,
                                        w_ym, w_yp, w_zm,
                                        w_ym_zm, w_yp_zm,
                                        dy);
        Grad_Vel[1][2] = dw_dy;
        
        amrex::Real dv_dy;
        FDM_GradientY_ZFaceFromCenter_D(dv_dy,
                                        v_c,
                                        v_ym, v_yp, v_zm,
                                        v_ym_zm, v_yp_zm,
                                        dy);
        Grad_Vel[1][1] = dv_dy;
        
        amrex::Real dv_dz;
        FDM_GradientZ_ZFaceFromCenter_D(dv_dz,
                                        v_c,
                                        v_zm,
                                        dz);
        Grad_Vel[2][1] = dv_dz;
        
        amrex::Real dw_dz;
        FDM_GradientZ_ZFaceFromCenter_D(dw_dz,
                                        w_c,
                                        w_zm,
                                        dz);
        Grad_Vel[2][2] = dw_dz;
        //Temperature gradient
        amrex::Real dT_dz;
        FDM_GradientZ_ZFaceFromCenter_D(dT_dz,
                                        T_c,
                                        T_zm,
                                        dz);
        Grad_T = dT_dz;
        //Material property
        //mu
        Interpolation2Face_Linear_D(mu,
                                    mu_zm,
                                    mu_c);
        //mub
        Interpolation2Face_Linear_D(mub,
                                    mub_zm,
                                    mub_c);
        //Kappa
        Interpolation2Face_Linear_D(Kappa,
                                    Kappa_zm,
                                    Kappa_c);
    }
#endif
    
    FluxDiffusion(Flux,
                  dmax,
                  tau,
                  Vel,
                  Grad_Vel,
                  Grad_T,
                  mu,
                  mub,
                  Kappa,
                  ID_Dimension);
    
}

#endif
#endif