// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO5B_D_Compressible6Eq_IE_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_WENO5B_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// Bounded WENO5 for the N-phase six-equation model, following the N-phase five-equation WENO5B.
// The partial densities and volume fractions share nonlinear weights and are limited so that
// they stay non-negative and the volume fractions sum to one. The velocities, the phasic
// Alpha_k p_k and the mixture pressure are reconstructed with WENO5.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO5B_D(amrex::Array<amrex::Real,NSTATE>& Conservative_Face,
                                                   amrex::Array<amrex::Real,NSTATE> const& Conservative_mm,
                                                   amrex::Array<amrex::Real,NSTATE> const& Conservative_m,
                                                   amrex::Array<amrex::Real,NSTATE> const& Conservative_c,
                                                   amrex::Array<amrex::Real,NSTATE> const& Conservative_p,
                                                   amrex::Array<amrex::Real,NSTATE> const& Conservative_pp,
                                                   Parm const& EOS_Parameter,
                                                   amrex::Real const& dx,
                                                   int ID_Face,
                                                   int ID_Poly)
{
    int const Num_Point   = 5;
    int const Num_Stencil = 3;
    int const i_c         = 2; //Index of the center cell in the stencil

    //Primitive variables on the stencil
    amrex::Array<amrex::Array<amrex::Real,NSTATE>,Num_Point> Conservative_s = {Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp};
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,NPHASE>         Alpha_rho_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,NPHASE>         Alpha_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,NPHASE>         Alpha_Ps_s;
    amrex::Array<amrex::Array<amrex::Real,Num_Point>,AMREX_SPACEDIM> Vel_s;
    amrex::Array<amrex::Real,Num_Point>                              p_s;
    for (int i_s = 0; i_s < Num_Point; i_s++){
        amrex::Array<amrex::Real,NSTATE> Primitive;
        Conservative2Primitive(Primitive,
                               Conservative_s[i_s],
                               EOS_Parameter);
        amrex::Array<amrex::Real,NPHASE> Alpha;
        amrex::Array<amrex::Real,NPHASE> Alpha_rho;
        amrex::Array<amrex::Real,NPHASE> Alpha_Ps;
        AMREX_D_TERM(amrex::Real u;,
                     amrex::Real v;,
                     amrex::Real w);
        amrex::Real p;
        Primitive2State(Alpha, Alpha_rho, Alpha_Ps,
           AMREX_D_DECL(u, v, w),
                        p,
                        Primitive);
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Alpha_rho_s[i_Phase][i_s] = Alpha_rho[i_Phase];
            Alpha_s    [i_Phase][i_s] = Alpha    [i_Phase];
            Alpha_Ps_s [i_Phase][i_s] = Alpha_Ps [i_Phase];
        }
        AMREX_D_TERM(Vel_s[0][i_s] = u;,
                     Vel_s[1][i_s] = v;,
                     Vel_s[2][i_s] = w;)
        p_s[i_s] = p;
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    amrex::Real Weight = (ID_Poly == 2) ? 1.0/6.0 : 1.0/12.0;
    int ID_BoundType = -1;

    //Mass and volume fraction, reconstructed at both faces of the center cell
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_rho_Weight;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Alpha_Weight;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoL_fp;
    amrex::Array<amrex::Real,NPHASE> AlphaL_fp;
    amrex::Array<amrex::Real,NPHASE> Alpha_rhoR_fm;
    amrex::Array<amrex::Real,NPHASE> AlphaR_fm;
    amrex::Real Sum_AlphaL_fp = 0.0;
    amrex::Real Sum_AlphaR_fm = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Array<amrex::Real,Num_Point> const& f_rho = Alpha_rho_s[i_Phase];
        amrex::Array<amrex::Real,Num_Point> const& f_a   = Alpha_s    [i_Phase];
        for (int i_Face = 0; i_Face < 2; i_Face++){
            int ID_FaceRe = (i_Face == 0) ? 1 : -1;
            amrex::Real Alpha_rho_Face;
            amrex::Real Alpha_Face;
            WENO5(Alpha_rho_WENO, Alpha_rho_Weight,
                  f_rho[0], f_rho[1], f_rho[2], f_rho[3], f_rho[4],
                  ID_FaceRe);
            Reconstruction2Face_WENO(Alpha_rho_Face,
                                     Alpha_rho_WENO, Alpha_rho_Weight,
                                     Num_Stencil);
            //The volume fraction uses the weights of the partial density
            WENO5(Alpha_WENO, Alpha_Weight,
                  f_a[0], f_a[1], f_a[2], f_a[3], f_a[4],
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
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
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

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_Face;
    amrex::Array<amrex::Real,NPHASE> Alpha_Face;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        f_LimiterBound(Alpha_rho_Face[i_Phase],
                       (ID_Face == 1) ? Alpha_rhoL_fp[i_Phase] : Alpha_rhoR_fm[i_Phase],
                       Alpha_rho_s[i_Phase][i_c],
                       Theta);
        f_LimiterBound(Alpha_Face[i_Phase],
                       (ID_Face == 1) ? AlphaL_fp[i_Phase] : AlphaR_fm[i_Phase],
                       Alpha_s[i_Phase][i_c],
                       Theta);
    }

    //Velocity, phasic pressures and mixture pressure
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel_Face;
    for (int i_Dim = 0; i_Dim < AMREX_SPACEDIM; i_Dim++){
        amrex::Array<amrex::Real,Num_Point> const& f = Vel_s[i_Dim];
        Reconstruction2Face_WENO5_D(Vel_Face[i_Dim],
                                    f[0], f[1], f[2], f[3], f[4],
                                    ID_Face);
    }
    amrex::Array<amrex::Real,NPHASE> Alpha_Ps_Face;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Array<amrex::Real,Num_Point> const& f = Alpha_Ps_s[i_Phase];
        Reconstruction2Face_WENO5_D(Alpha_Ps_Face[i_Phase],
                                    f[0], f[1], f[2], f[3], f[4],
                                    ID_Face);
    }
    amrex::Real p_Face;
    Reconstruction2Face_WENO5_D(p_Face,
                                p_s[0], p_s[1], p_s[2], p_s[3], p_s[4],
                                ID_Face);

    amrex::Array<amrex::Real,NSTATE> Primitive_Face;
    State2Primitive(Primitive_Face,
                    Alpha_Face, Alpha_rho_Face, Alpha_Ps_Face,
       AMREX_D_DECL(Vel_Face[0], Vel_Face[1], Vel_Face[2]),
                    p_Face);
    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
}

#endif
#endif
