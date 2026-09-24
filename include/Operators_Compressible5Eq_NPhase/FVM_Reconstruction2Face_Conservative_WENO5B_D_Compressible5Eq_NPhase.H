// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO5B_D_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_WENO5B_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO5B_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Conservative_Face,
                                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_mm,
                                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_m,
                                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_p,
                                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_pp,
                                                   Parm const& EOS_Parameter,
                                                   amrex::Real const& dx,
                                                   int ID_Face,
                                                   int ID_Poly)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_mm;
    Conservative2Primitive(Primitive_mm,
                           Conservative_mm,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_mm;
    AMREX_D_TERM(amrex::Real u_mm;,
                 amrex::Real v_mm;,
                 amrex::Real w_mm);
    amrex::Real p_mm;
    amrex::Array<amrex::Real,NPHASE> Alpha_mm;
    Primitive2State(Alpha_rho_mm,
       AMREX_D_DECL(u_mm, v_mm, w_mm),
                    p_mm,
                    Alpha_mm,
                    Primitive_mm);

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_m;
    Conservative2Primitive(Primitive_m,
                           Conservative_m,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_m;
    AMREX_D_TERM(amrex::Real u_m;,
                 amrex::Real v_m;,
                 amrex::Real w_m);
    amrex::Real p_m;
    amrex::Array<amrex::Real,NPHASE> Alpha_m;
    Primitive2State(Alpha_rho_m,
       AMREX_D_DECL(u_m, v_m, w_m),
                    p_m,
                    Alpha_m,
                    Primitive_m);

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_c;
    AMREX_D_TERM(amrex::Real u_c;,
                 amrex::Real v_c;,
                 amrex::Real w_c);
    amrex::Real p_c;
    amrex::Array<amrex::Real,NPHASE> Alpha_c;
    Primitive2State(Alpha_rho_c,
       AMREX_D_DECL(u_c, v_c, w_c),
                    p_c,
                    Alpha_c,
                    Primitive_c);

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_p;
    Conservative2Primitive(Primitive_p,
                           Conservative_p,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_p;
    AMREX_D_TERM(amrex::Real u_p;,
                 amrex::Real v_p;,
                 amrex::Real w_p);
    amrex::Real p_p;
    amrex::Array<amrex::Real,NPHASE> Alpha_p;
    Primitive2State(Alpha_rho_p,
       AMREX_D_DECL(u_p, v_p, w_p),
                    p_p,
                    Alpha_p,
                    Primitive_p);

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_pp;
    Conservative2Primitive(Primitive_pp,
                           Conservative_pp,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_pp;
    AMREX_D_TERM(amrex::Real u_pp;,
                 amrex::Real v_pp;,
                 amrex::Real w_pp);
    amrex::Real p_pp;
    amrex::Array<amrex::Real,NPHASE> Alpha_pp;
    Primitive2State(Alpha_rho_pp,
       AMREX_D_DECL(u_pp, v_pp, w_pp),
                    p_pp,
                    Alpha_pp,
                    Primitive_pp);

    //Mass and Volume Fraction
    int const Num_Stencil = 3;
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    amrex::Real Weight;
    if (ID_Poly == 2){
        Weight = 1.0/6.0;
    }
    else if (ID_Poly == 4){
        Weight = 1.0/12.0;
    }

    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_Weight;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_Weight;

    int ID_FaceRe;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL_fp;
    amrex::Array<amrex::Real,NPHASE> AlphaL_fp;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR_fm;
    amrex::Array<amrex::Real,NPHASE> AlphaR_fm;
    amrex::Real Sum_AlphaL_fp = 0.0;
    amrex::Real Sum_AlphaR_fm = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        //fL_fp
        ID_FaceRe = 1;
        //Alpha_rho
        WENO5(Alpha_rho_WENO, Alpha_rho_Weight,
              Alpha_rho_mm[i_Phase], Alpha_rho_m[i_Phase], Alpha_rho_c[i_Phase], Alpha_rho_p[i_Phase], Alpha_rho_pp[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(Alpha_rhoL_fp[i_Phase],
                                 Alpha_rho_WENO, Alpha_rho_Weight,
                                 Num_Stencil);
        //Alpha
        WENO5(Alpha_WENO, Alpha_Weight,
              Alpha_mm[i_Phase], Alpha_m[i_Phase], Alpha_c[i_Phase], Alpha_p[i_Phase], Alpha_pp[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(AlphaL_fp[i_Phase],
                                 Alpha_WENO, Alpha_rho_Weight,
                                 Num_Stencil);

        //fR_fm
        ID_FaceRe = -1;
        //Alpha_rho
        WENO5(Alpha_rho_WENO, Alpha_rho_Weight,
              Alpha_rho_mm[i_Phase], Alpha_rho_m[i_Phase], Alpha_rho_c[i_Phase], Alpha_rho_p[i_Phase], Alpha_rho_pp[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(Alpha_rhoR_fm[i_Phase],
                                 Alpha_rho_WENO, Alpha_rho_Weight,
                                 Num_Stencil);
        //Alpha
        WENO5(Alpha_WENO, Alpha_Weight,
              Alpha_mm[i_Phase], Alpha_m[i_Phase], Alpha_c[i_Phase], Alpha_p[i_Phase], Alpha_pp[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(AlphaR_fm[i_Phase],
                                 Alpha_WENO, Alpha_rho_Weight,
                                 Num_Stencil);

        //Positivity
        int ID_BoundType = -1;
        amrex::Real Theta_Alpha_rho;
        Limiter_Bound(Theta_Alpha_rho,
                      inf, 0.0,
                      Alpha_rho_c[i_Phase], Alpha_rhoL_fp[i_Phase], Alpha_rhoR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta_Alpha;
        Limiter_Bound(Theta_Alpha,
                      inf, 0.0,
                      Alpha_c[i_Phase], AlphaL_fp[i_Phase], AlphaR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta = std::min(Theta_Alpha_rho, Theta_Alpha);

        amrex::Real f_Limit;
        //Alpha_rhoL_fp
        f_LimiterBound(f_Limit,
                       Alpha_rhoL_fp[i_Phase], Alpha_rho_c[i_Phase],
                       Theta);
        Alpha_rhoL_fp[i_Phase] = f_Limit;
        //AlphaL_fp
        f_LimiterBound(f_Limit,
                       AlphaL_fp[i_Phase], Alpha_c[i_Phase],
                       Theta);
        AlphaL_fp[i_Phase] = f_Limit;
        //Alpha_rhoR_fm
        f_LimiterBound(f_Limit,
                       Alpha_rhoR_fm[i_Phase], Alpha_rho_c[i_Phase],
                       Theta);
        Alpha_rhoR_fm[i_Phase] = f_Limit;
        //AlphaR_fm
        f_LimiterBound(f_Limit,
                       AlphaR_fm[i_Phase], Alpha_c[i_Phase],
                       Theta);
        AlphaR_fm[i_Phase] = f_Limit;

        //Summation
        Sum_AlphaL_fp += AlphaL_fp[i_Phase];
        Sum_AlphaR_fm += AlphaR_fm[i_Phase];
    }
    //Summation and Positivity
    amrex::Real Theta = 1.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        //Summation
        Alpha_rhoL_fp[i_Phase] /= Sum_AlphaL_fp;
        AlphaL_fp    [i_Phase] /= Sum_AlphaL_fp;

        Alpha_rhoR_fm[i_Phase] /= Sum_AlphaR_fm;
        AlphaR_fm    [i_Phase] /= Sum_AlphaR_fm;

        //Positivity
        int ID_BoundType = -1;
        amrex::Real Theta_Alpha_rho;
        Limiter_Bound(Theta_Alpha_rho,
                      inf, 0.0,
                      Alpha_rho_c[i_Phase], Alpha_rhoL_fp[i_Phase], Alpha_rhoR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta_Alpha;
        Limiter_Bound(Theta_Alpha,
                      inf, 0.0,
                      Alpha_c[i_Phase], AlphaL_fp[i_Phase], AlphaR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        Theta = std::min({Theta, Theta_Alpha_rho, Theta_Alpha});
    }
    //Alpha_rho_Face and Alpha_Face
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Face;
    amrex::Array<amrex::Real,NPHASE> Alpha_Face;
    //fL_fp
    if (ID_Face == 1){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            //Alpha_rho
            f_LimiterBound(Alpha_rho_Face[i_Phase],
                           Alpha_rhoL_fp[i_Phase], Alpha_rho_c[i_Phase],
                           Theta);
            //Alpha
            f_LimiterBound(Alpha_Face[i_Phase],
                           AlphaL_fp[i_Phase], Alpha_c[i_Phase],
                           Theta);
        }
    }
    //fR_fm
    else if (ID_Face == -1){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            //Alpha_rho
            f_LimiterBound(Alpha_rho_Face[i_Phase],
                           Alpha_rhoR_fm[i_Phase], Alpha_rho_c[i_Phase],
                           Theta);
            //Alpha
            f_LimiterBound(Alpha_Face[i_Phase],
                           AlphaR_fm[i_Phase], Alpha_c[i_Phase],
                           Theta);
        }
    }

    //Velocity
    amrex::Real u_Face;
    Reconstruction2Face_WENO5_D(u_Face,
                                u_mm, u_m, u_c, u_p, u_pp,
                                ID_Face);
#if (AMREX_SPACEDIM > 1)
    amrex::Real v_Face;
    Reconstruction2Face_WENO5_D(v_Face,
                                v_mm, v_m, v_c, v_p, v_pp,
                                ID_Face);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real w_Face;
    Reconstruction2Face_WENO5_D(w_Face,
                                w_mm, w_m, w_c, w_p, w_pp,
                                ID_Face);
#endif
    //Pressure
    amrex::Real p_Face;
    Reconstruction2Face_WENO5_D(p_Face,
                                p_mm, p_m, p_c, p_p, p_pp,
                                ID_Face);

    //Conservative_Face
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive_Face;
    State2Primitive(Primitive_Face,
                    Alpha_rho_Face,
       AMREX_D_DECL(u_Face, v_Face, w_Face),
                    p_Face,
                    Alpha_Face);
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif