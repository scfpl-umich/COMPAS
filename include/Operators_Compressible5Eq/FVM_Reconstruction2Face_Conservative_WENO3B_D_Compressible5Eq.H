// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO3B_D_Compressible5Eq_H_
#define FVM_Reconstruction2Face_Conservative_WENO3B_D_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO3B_D(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Conservative_Face,
                                                  amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_m,
                                                  amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_c,
                                                  amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_p,
                                                  Parm const& EOS_Parameter,
                                                  int ID_Face)
{
    int const Num_Phase = 2;

    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_m;
    Conservative2Primitive(Primitive_m,
                           Conservative_m,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1_m;
    amrex::Real Alpha2_rho2_m;
    AMREX_D_TERM(amrex::Real u_m;,
                 amrex::Real v_m;,
                 amrex::Real w_m);
    amrex::Real p_m;
    amrex::Real Alpha1_m;
    Primitive2State(Alpha1_rho1_m, Alpha2_rho2_m,
       AMREX_D_DECL(u_m, v_m, w_m),
                    p_m,
                    Alpha1_m,
                    Primitive_m);
    amrex::Real Alpha2_m = 1.0 - Alpha1_m;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rho_m = {Alpha1_rho1_m, Alpha2_rho2_m};
    amrex::Array<amrex::Real,Num_Phase> Alpha_m     = {Alpha1_m     , Alpha2_m     };
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1_c;
    amrex::Real Alpha2_rho2_c;
    AMREX_D_TERM(amrex::Real u_c;,
                 amrex::Real v_c;,
                 amrex::Real w_c);
    amrex::Real p_c;
    amrex::Real Alpha1_c;
    Primitive2State(Alpha1_rho1_c, Alpha2_rho2_c,
       AMREX_D_DECL(u_c, v_c, w_c),
                    p_c,
                    Alpha1_c,
                    Primitive_c);
    amrex::Real Alpha2_c = 1.0 - Alpha1_c;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rho_c = {Alpha1_rho1_c, Alpha2_rho2_c};
    amrex::Array<amrex::Real,Num_Phase> Alpha_c     = {Alpha1_c     , Alpha2_c     };
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_p;
    Conservative2Primitive(Primitive_p,
                           Conservative_p,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1_p;
    amrex::Real Alpha2_rho2_p;
    AMREX_D_TERM(amrex::Real u_p;,
                 amrex::Real v_p;,
                 amrex::Real w_p);
    amrex::Real p_p;
    amrex::Real Alpha1_p;
    Primitive2State(Alpha1_rho1_p, Alpha2_rho2_p,
       AMREX_D_DECL(u_p, v_p, w_p),
                    p_p,
                    Alpha1_p,
                    Primitive_p);
    amrex::Real Alpha2_p = 1.0 - Alpha1_p;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rho_p = {Alpha1_rho1_p, Alpha2_rho2_p};
    amrex::Array<amrex::Real,Num_Phase> Alpha_p     = {Alpha1_p     , Alpha2_p     };
    
    //Mass and Volume Fraction
    int const Num_Stencil = 2;
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    amrex::Real Weight = 1.0/6.0;

    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_Weight;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_Weight;

    int ID_FaceRe;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rhoL_fp;
    amrex::Array<amrex::Real,Num_Phase> AlphaL_fp;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rhoR_fm;
    amrex::Array<amrex::Real,Num_Phase> AlphaR_fm;
    amrex::Real Sum_AlphaL_fp = 0.0;
    amrex::Real Sum_AlphaR_fm = 0.0;
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        //fL_fp
        ID_FaceRe = 1;
        //Alpha_rho
        WENO3(Alpha_rho_WENO, Alpha_rho_Weight,
              Alpha_rho_m[i_Phase], Alpha_rho_c[i_Phase], Alpha_rho_p[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(Alpha_rhoL_fp[i_Phase],
                                 Alpha_rho_WENO, Alpha_rho_Weight,
                                 Num_Stencil);
        //Alpha
        WENO3(Alpha_WENO, Alpha_Weight,
              Alpha_m[i_Phase], Alpha_c[i_Phase], Alpha_p[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(AlphaL_fp[i_Phase],
                                 Alpha_WENO, Alpha_rho_Weight,
                                 Num_Stencil);

        //fR_fm
        ID_FaceRe = -1;
        //Alpha_rho
        WENO3(Alpha_rho_WENO, Alpha_rho_Weight,
              Alpha_rho_m[i_Phase], Alpha_rho_c[i_Phase], Alpha_rho_p[i_Phase],
              ID_FaceRe);
        Reconstruction2Face_WENO(Alpha_rhoR_fm[i_Phase],
                                 Alpha_rho_WENO, Alpha_rho_Weight,
                                 Num_Stencil);
        //Alpha
        WENO3(Alpha_WENO, Alpha_Weight,
              Alpha_m[i_Phase], Alpha_c[i_Phase], Alpha_p[i_Phase],
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
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
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
    amrex::Real Alpha1_rho1_Face;
    amrex::Real Alpha2_rho2_Face;
    amrex::Real Alpha1_Face;
    //fL_fp
    if (ID_Face == 1){
        Alpha1_rho1_Face = Alpha_rhoL_fp[0];
        Alpha2_rho2_Face = Alpha_rhoL_fp[1];
        Alpha1_Face      = AlphaL_fp    [0];
    }
    //fR_fm
    else if (ID_Face == -1){
        Alpha1_rho1_Face = Alpha_rhoR_fm[0];
        Alpha2_rho2_Face = Alpha_rhoR_fm[1];
        Alpha1_Face      = AlphaR_fm    [0];
    }
    amrex::Real f_Limit;
    //Alpha1_rho1_Face
    f_LimiterBound(f_Limit,
                   Alpha1_rho1_Face, Alpha1_rho1_c,
                   Theta);
    Alpha1_rho1_Face = f_Limit;
    //Alpha2_rho2_Face
    f_LimiterBound(f_Limit,
                   Alpha2_rho2_Face, Alpha2_rho2_c,
                   Theta);
    Alpha2_rho2_Face = f_Limit;
    //Alpha1_Face
    f_LimiterBound(f_Limit,
                   Alpha1_Face, Alpha1_c,
                   Theta);
    Alpha1_Face = f_Limit;
    
    //Velocity
    amrex::Real u_Face;
    Reconstruction2Face_WENO3_D(u_Face,
                                u_m, u_c, u_p,
                                ID_Face);
#if (AMREX_SPACEDIM > 1)
    amrex::Real v_Face;
    Reconstruction2Face_WENO3_D(v_Face,
                                v_m, v_c, v_p,
                                ID_Face);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real w_Face;
    Reconstruction2Face_WENO3_D(w_Face,
                                w_m, w_c, w_p,
                                ID_Face);
#endif
    //Pressure
    amrex::Real p_Face;
    Reconstruction2Face_WENO3_D(p_Face,
                                p_m, p_c, p_p,
                                ID_Face);

    //Conservative_Face
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive_Face;
    State2Primitive(Primitive_Face,
                    Alpha1_rho1_Face, Alpha2_rho2_Face,
       AMREX_D_DECL(u_Face, v_Face, w_Face),
                    p_Face,
                    Alpha1_Face);
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif