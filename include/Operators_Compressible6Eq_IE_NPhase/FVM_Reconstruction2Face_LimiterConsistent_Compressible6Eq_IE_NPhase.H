// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_LimiterConsistent_Compressible6Eq_IE_NPhase_H_
#define FVM_Reconstruction2Face_LimiterConsistent_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_LimiterConsistent(amrex::Array<amrex::Real,NPHASE>& Alpha_Ps_Limit,
                                               amrex::Array<amrex::Real,NPHASE>& Alpha_rho_Limit,
                                               amrex::Array<amrex::Real,NPHASE>& Alpha_Limit,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_Ps_Face,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_Face,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_Face,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_Ps_c,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_c,
                                               amrex::Array<amrex::Real,NPHASE> const& Alpha_c,
                                               amrex::Real const& dx,
                                               int ID_Face)
{
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha_Ps;
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Slope_Alpha;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real f_c;
        amrex::Real f_Face;
        //Alpha_Ps
        f_c    = Alpha_Ps_c   [i_Phase];
        f_Face = Alpha_Ps_Face[i_Phase];
        Slope_Alpha_Ps[i_Phase] = (f_Face - f_c)/(amrex::Real(ID_Face)*0.5*dx);

        //Alpha_rho
        f_c    = Alpha_rho_c   [i_Phase];
        f_Face = Alpha_rho_Face[i_Phase];
        Slope_Alpha_rho[i_Phase] = (f_Face - f_c)/(amrex::Real(ID_Face)*0.5*dx);
        
        //Alpha
        f_c    = Alpha_c   [i_Phase];
        f_Face = Alpha_Face[i_Phase];
        Slope_Alpha[i_Phase] = (f_Face - f_c)/(amrex::Real(ID_Face)*0.5*dx);
    }
    amrex::Array<amrex::Real,NPHASE> Theta;
    Limiter_Consistent(Theta,
                       Slope_Alpha);

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real f_c;
        amrex::Real Slope;
        amrex::Real f_Face;
        //Alpha_P
        f_c = Alpha_Ps_c[i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha_Ps[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   Slope,
                                   dx,
                                   ID_Face);
        Alpha_Ps_Limit[i_Phase] = f_Face;

        //Alpha_rho
        f_c = Alpha_rho_c[i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha_rho[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   Slope,
                                   dx,
                                   ID_Face);
        Alpha_rho_Limit[i_Phase] = f_Face;
        
        //Alpha
        f_c = Alpha_c[i_Phase];
        Slope = Theta[i_Phase]*Slope_Alpha[i_Phase];
        Reconstruction2Face_Linear(f_Face,
                                   f_c,
                                   Slope,
                                   dx,
                                   ID_Face);
        Alpha_Limit[i_Phase] = f_Face;
    }
    
}

#endif
#endif