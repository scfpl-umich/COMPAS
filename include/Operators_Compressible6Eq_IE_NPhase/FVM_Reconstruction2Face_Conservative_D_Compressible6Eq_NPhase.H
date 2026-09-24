// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_D_Compressible6Eq_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_D_Compressible6Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_D(amrex::Array<amrex::Real,NSTATE>& Conservative_Face,
                                            amrex::Array<amrex::Real,NSTATE> const& Conservative_mm,
                                            amrex::Array<amrex::Real,NSTATE> const& Conservative_m,
                                            amrex::Array<amrex::Real,NSTATE> const& Conservative_c,
                                            amrex::Array<amrex::Real,NSTATE> const& Conservative_p,
                                            amrex::Array<amrex::Real,NSTATE> const& Conservative_pp,
                                            Parm const& EOS_Parameter,
                                            amrex::Real const& dx,
                                            int ID_Face,
                                            amrex::Array<int,NSTATE> const& ID_WENO,
                                            amrex::Array<int,NSTATE> const& ID_Limiter,
                                            amrex::Array<int,NSTATE> const& ID_Poly)
{
    amrex::Array<amrex::Real,NSTATE> Primitive_mm;
    Conservative2Primitive(Primitive_mm,
                           Conservative_mm,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_m;
    Conservative2Primitive(Primitive_m,
                           Conservative_m,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_c;
    Conservative2Primitive(Primitive_c,
                           Conservative_c,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_p;
    Conservative2Primitive(Primitive_p,
                           Conservative_p,
                           EOS_Parameter);
    amrex::Array<amrex::Real,NSTATE> Primitive_pp;
    Conservative2Primitive(Primitive_pp,
                           Conservative_pp,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,NSTATE> Primitive_Face;
    for (int i_State = 0; i_State < NSTATE; i_State++){
        amrex::Real f_mm = Primitive_mm[i_State];
        amrex::Real f_m  = Primitive_m [i_State];
        amrex::Real f_c  = Primitive_c [i_State];
        amrex::Real f_p  = Primitive_p [i_State];
        amrex::Real f_pp = Primitive_pp[i_State];
        
        amrex::Real f_Face;
        amrex::Real f_Stencil[] = {f_mm, f_m, f_c, f_p, f_pp};
        Reconstruction2Face_D(f_Face,
                              f_Stencil,
                              EOS_Parameter.UpperBound[i_State],
                              EOS_Parameter.LowerBound[i_State],
                              dx,
                              ID_Face,
                              ID_WENO[i_State],
                              ID_Limiter[i_State],
                              ID_Poly[i_State],
                              EOS_Parameter.ID_BoundType[i_State]);
        Primitive_Face[i_State] = f_Face;
    }


    // consistent alpha
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha;
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha_rho_e;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real f_c;
        amrex::Real f_Face;

        //Alpha
        f_c = Primitive_c[i_Phase];
        f_Face = Primitive_Face[i_Phase];
        Slope_Alpha[i_Phase] = (f_Face - f_c)/(amrex::Real(ID_Face)*0.5*dx);


        //Alpha_rho
        f_c = Primitive_c[NPHASE + i_Phase];
        f_Face = Primitive_Face[NPHASE + i_Phase];
        Slope_Alpha_rho[i_Phase] = (f_Face -f_c)/(amrex::Real(ID_Face)*0.5*dx);

        //Alpha_rho_e
        f_c = Primitive_c[2*NPHASE + i_Phase];
        f_Face = Primitive_Face[2*NPHASE + i_Phase];
        Slope_Alpha_rho_e[i_Phase] = (f_Face -f_c)/(amrex::Real(ID_Face)*0.5*dx);

    }

    amrex::Array<amrex::Real,NPHASE> Theta;
    Limiter_Consistent(Theta, Slope_Alpha);

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real f_c;
        amrex::Real Slope;
        amrex::Real f_Face;

        //Alpha
        f_c = Primitive_c[i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   dx,
                                   Slope,
                                   ID_Face);
        Primitive_Face[i_Phase] = f_Face;

        //Alpha_rho
        f_c = Primitive_c[NPHASE+i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha_rho[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   dx,
                                   Slope,
                                   ID_Face);
        Primitive_Face[NPHASE+i_Phase] = f_Face;

        //Alpha_rho_e
        f_c = Primitive_c[ 2*NPHASE + i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha_rho_e[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   dx,
                                   Slope,
                                   ID_Face);
        Primitive_Face[ 2*NPHASE +i_Phase] = f_Face;
        
    }
    //


    Primitive2Conservative(Conservative_Face,
                           Primitive_Face,
                           EOS_Parameter);
    
}

#endif
#endif