// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO3B_D_Compressible6Eq_H_
#define FVM_Reconstruction2Face_Conservative_WENO3B_D_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

// Bounded WENO3 for the two-phase six-equation model, following the five-equation WENO3B.
// The partial densities and volume fractions share nonlinear weights and are limited so that
// they stay non-negative and the volume fractions sum to one. The velocities and the phasic
// Alpha_k p_k are reconstructed with WENO3.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO3B_D(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Conservative_Face,
                                                   amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative_m,
                                                   amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative_c,
                                                   amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative_p,
                                                   Parm const& EOS_Parameter,
                                                   amrex::Real const& dx,
                                                   int ID_Face)
{
    int const Num_Phase   = 2;
    int const Num_Point   = 3;
    int const Num_Stencil = 2;
    int const i_c         = 1; //Index of the center cell in the stencil

    //Primitive variables on the stencil
    amrex::Array<amrex::Array<amrex::Real,5+AMREX_SPACEDIM>,Num_Point> Conservative_s = {Conservative_m, Conservative_c, Conservative_p};
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,Num_Phase>      Alpha_rho_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,Num_Phase>      Alpha_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,Num_Phase>      Alpha_p_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,AMREX_SPACEDIM> Vel_s;
    for (int i_s = 0; i_s < Num_Point; i_s++){
        amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive;
        Conservative2Primitive(Primitive,
                               Conservative_s[i_s],
                               EOS_Parameter);
        amrex::Real Alpha1_rho1;
        amrex::Real Alpha2_rho2;
        AMREX_D_TERM(amrex::Real u;,
                     amrex::Real v;,
                     amrex::Real w);
        amrex::Real Alpha1_p1;
        amrex::Real Alpha2_p2;
        amrex::Real Alpha1;
        Primitive2State(Alpha1_rho1, Alpha2_rho2,
           AMREX_D_DECL(u, v, w),
                        Alpha1_p1, Alpha2_p2,
                        Alpha1,
                        Primitive);
        Alpha_rho_s[0][i_s] = Alpha1_rho1;
        Alpha_rho_s[1][i_s] = Alpha2_rho2;
        Alpha_s    [0][i_s] = Alpha1;
        Alpha_s    [1][i_s] = 1.0 - Alpha1;
        Alpha_p_s  [0][i_s] = Alpha1_p1;
        Alpha_p_s  [1][i_s] = Alpha2_p2;
        AMREX_D_TERM(Vel_s[0][i_s] = u;,
                     Vel_s[1][i_s] = v;,
                     Vel_s[2][i_s] = w;)
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    amrex::Real Weight = 1.0/6.0;
    int ID_BoundType = -1;

    //Mass and volume fraction, reconstructed at both faces of the center cell
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_Weight;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_Weight;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rhoL_fp;
    amrex::Array<amrex::Real,Num_Phase> AlphaL_fp;
    amrex::Array<amrex::Real,Num_Phase> Alpha_rhoR_fm;
    amrex::Array<amrex::Real,Num_Phase> AlphaR_fm;
    amrex::Real Sum_AlphaL_fp = 0.0;
    amrex::Real Sum_AlphaR_fm = 0.0;
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        amrex::Array<amrex::Real,Num_Point> const& f_rho = Alpha_rho_s[i_Phase];
        amrex::Array<amrex::Real,Num_Point> const& f_a   = Alpha_s    [i_Phase];
        for (int i_Face = 0; i_Face < 2; i_Face++){
            int ID_FaceRe = (i_Face == 0) ? 1 : -1;
            amrex::Real Alpha_rho_Face;
            amrex::Real Alpha_Face;
            WENO3(Alpha_rho_WENO, Alpha_rho_Weight,
                  f_rho[0], f_rho[1], f_rho[2],
                  ID_FaceRe);
            Reconstruction2Face_WENO(Alpha_rho_Face,
                                     Alpha_rho_WENO, Alpha_rho_Weight,
                                     Num_Stencil);
            //The volume fraction uses the weights of the partial density
            WENO3(Alpha_WENO, Alpha_Weight,
                  f_a[0], f_a[1], f_a[2],
                  ID_FaceRe);
            Reconstruction2Face_WENO(Alpha_Face,
                                     Alpha_WENO, Alpha_rho_Weight,
                                     Num_Stencil);
            if (ID_FaceRe == 1){
                Alpha_rhoL_fp[i_Phase] = Alpha_rho_Face;
                AlphaL_fp    [i_Phase] = Alpha_Face;
            }
            else {
                Alpha_rhoR_fm[i_Phase] = Alpha_rho_Face;
                AlphaR_fm    [i_Phase] = Alpha_Face;
            }
        }

        //Positivity
        amrex::Real Theta_Alpha_rho;
        Limiter_Bound(Theta_Alpha_rho,
                      inf, 0.0,
                      f_rho[i_c], Alpha_rhoL_fp[i_Phase], Alpha_rhoR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta_Alpha;
        Limiter_Bound(Theta_Alpha,
                      inf, 0.0,
                      f_a[i_c], AlphaL_fp[i_Phase], AlphaR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta = std::min(Theta_Alpha_rho, Theta_Alpha);
        amrex::Real f_Limit;
        f_LimiterBound(f_Limit, Alpha_rhoL_fp[i_Phase], f_rho[i_c], Theta);
        Alpha_rhoL_fp[i_Phase] = f_Limit;
        f_LimiterBound(f_Limit, AlphaL_fp[i_Phase], f_a[i_c], Theta);
        AlphaL_fp[i_Phase] = f_Limit;
        f_LimiterBound(f_Limit, Alpha_rhoR_fm[i_Phase], f_rho[i_c], Theta);
        Alpha_rhoR_fm[i_Phase] = f_Limit;
        f_LimiterBound(f_Limit, AlphaR_fm[i_Phase], f_a[i_c], Theta);
        AlphaR_fm[i_Phase] = f_Limit;

        Sum_AlphaL_fp += AlphaL_fp[i_Phase];
        Sum_AlphaR_fm += AlphaR_fm[i_Phase];
    }

    //Summation: rescale so that the volume fractions sum to one, then limit again
    amrex::Real Theta = 1.0;
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        Alpha_rhoL_fp[i_Phase] /= Sum_AlphaL_fp;
        AlphaL_fp    [i_Phase] /= Sum_AlphaL_fp;
        Alpha_rhoR_fm[i_Phase] /= Sum_AlphaR_fm;
        AlphaR_fm    [i_Phase] /= Sum_AlphaR_fm;
        amrex::Real Theta_Alpha_rho;
        Limiter_Bound(Theta_Alpha_rho,
                      inf, 0.0,
                      Alpha_rho_s[i_Phase][i_c], Alpha_rhoL_fp[i_Phase], Alpha_rhoR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        amrex::Real Theta_Alpha;
        Limiter_Bound(Theta_Alpha,
                      inf, 0.0,
                      Alpha_s[i_Phase][i_c], AlphaL_fp[i_Phase], AlphaR_fm[i_Phase],
                      Weight,
                      ID_BoundType);
        Theta = std::min({Theta, Theta_Alpha_rho, Theta_Alpha});
    }

    amrex::Real Alpha1_rho1_Face;
    amrex::Real Alpha2_rho2_Face;
    amrex::Real Alpha1_Face;
    if (ID_Face == 1){
        Alpha1_rho1_Face = Alpha_rhoL_fp[0];
        Alpha2_rho2_Face = Alpha_rhoL_fp[1];
        Alpha1_Face      = AlphaL_fp    [0];
    }
    else {
        Alpha1_rho1_Face = Alpha_rhoR_fm[0];
        Alpha2_rho2_Face = Alpha_rhoR_fm[1];
        Alpha1_Face      = AlphaR_fm    [0];
    }
    amrex::Real f_Limit;
    f_LimiterBound(f_Limit, Alpha1_rho1_Face, Alpha_rho_s[0][i_c], Theta);
    Alpha1_rho1_Face = f_Limit;
    f_LimiterBound(f_Limit, Alpha2_rho2_Face, Alpha_rho_s[1][i_c], Theta);
    Alpha2_rho2_Face = f_Limit;
    f_LimiterBound(f_Limit, Alpha1_Face, Alpha_s[0][i_c], Theta);
    Alpha1_Face = f_Limit;

    //Velocity and phasic pressures
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel_Face;
    for (int i_Dim = 0; i_Dim < AMREX_SPACEDIM; i_Dim++){
        amrex::Array<amrex::Real,Num_Point> const& f = Vel_s[i_Dim];
        Reconstruction2Face_WENO3_D(Vel_Face[i_Dim],
                                    f[0], f[1], f[2],
                                    ID_Face);
    }
    amrex::Array<amrex::Real,Num_Phase> Alpha_p_Face;
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        amrex::Array<amrex::Real,Num_Point> const& f = Alpha_p_s[i_Phase];
        Reconstruction2Face_WENO3_D(Alpha_p_Face[i_Phase],
                                    f[0], f[1], f[2],
                                    ID_Face);
    }

    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive_Face;
    State2Primitive(Primitive_Face,
                    Alpha1_rho1_Face, Alpha2_rho2_Face,
       AMREX_D_DECL(Vel_Face[0], Vel_Face[1], Vel_Face[2]),
                    Alpha_p_Face[0], Alpha_p_Face[1],
                    Alpha1_Face);
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
}

#endif
#endif
