// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SurfaceTension_Pressure_K_H_
#define SurfaceTension_Pressure_K_H_

#include<Macros.H>
#if (SURFACE_TENSION == true)

// =======================================
// Mixture pressure of a cell and its replacement, for the initial pressure projection
// (Compressible_PhaseField::SurfaceTension_InitialPressure).
// =======================================

// Mixture pressure: p of the five-equation models, alpha_1 p_1 + alpha_2 p_2 of SIXEQS, and
// P = sum_k alpha_k p_k of SIXEQS_IE_NPHASE
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real SurfaceTension_Pressure_K(int i, int j, int k,
                                      amrex::Array4<amrex::Real const> const& U,
                                      Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NSTATE> Conservative;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = U(i,j,k,iState);
    }
    amrex::Array<amrex::Real,NSTATE> Primitive;
    Conservative2Primitive(Primitive, Conservative, EOS_Parameter);
#if (PHYSICS == FIVEEQS || PHYSICS == FIVEEQS_NPHASE)
    return Primitive[INDEX_Energy];
#elif (PHYSICS == SIXEQS)
    return Primitive[INDEX_Energy1] + Primitive[INDEX_Energy2];
#else
    return Primitive[INDEX_Total_Energy];
#endif
}

// Set the pressure of cell (i,j,k) to p, with the volume fractions, partial densities and velocity
// fixed: only the energy components change (SIXEQS: alpha_k p_k = alpha_k p; SIXEQS_IE_NPHASE:
// alpha_k p_k = alpha_k p and P = p). Returns 1 if p + pinf_k <= 0 for a present phase
// (alpha_k > 0), with pinf_k = B_k/(1 + A_k), and then leaves the cell unchanged.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int SurfaceTension_SetPressure_K(int i, int j, int k,
                                 amrex::Array4<amrex::Real> const& U,
                                 amrex::Real p,
                                 Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NSTATE> Conservative;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = U(i,j,k,iState);
    }
    amrex::Array<amrex::Real,NSTATE> Primitive;
    Conservative2Primitive(Primitive, Conservative, EOS_Parameter);

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    amrex::Real Alpha1 = Conservative[INDEX_VolumeFraction1];
    amrex::Array<amrex::Real,2> Alpha = {Alpha1, 1.0 - Alpha1};
    amrex::Array<amrex::Real,2> pinf  = {EOS_Parameter.B1/(1.0 + EOS_Parameter.A1),
                                         EOS_Parameter.B2/(1.0 + EOS_Parameter.A2)};
    for (int i_Phase = 0; i_Phase < 2; i_Phase++){
        if (Alpha[i_Phase] > 0.0 && !(p + pinf[i_Phase] > 0.0)){ return 1; }
    }
#else
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real pinf = EOS_Parameter.B[i_Phase]/(1.0 + EOS_Parameter.A[i_Phase]);
        if (Conservative[INDEX_VolumeFraction1+i_Phase] > 0.0 && !(p + pinf > 0.0)){ return 1; }
    }
#endif

#if (PHYSICS == FIVEEQS || PHYSICS == FIVEEQS_NPHASE)
    Primitive[INDEX_Energy] = p;
#elif (PHYSICS == SIXEQS)
    Primitive[INDEX_Energy1] = Alpha[0]*p;
    Primitive[INDEX_Energy2] = Alpha[1]*p;
#else
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Primitive[INDEX_Energy1+i_Phase] = Conservative[INDEX_VolumeFraction1+i_Phase]*p;
    }
    Primitive[INDEX_Total_Energy] = p;
#endif

    amrex::Array<amrex::Real,NSTATE> Rebuilt;
    Primitive2Conservative(Rebuilt, Primitive, EOS_Parameter);

#if (PHYSICS == FIVEEQS || PHYSICS == FIVEEQS_NPHASE)
    U(i,j,k,INDEX_Energy) = Rebuilt[INDEX_Energy];
#elif (PHYSICS == SIXEQS)
    U(i,j,k,INDEX_Energy1) = Rebuilt[INDEX_Energy1];
    U(i,j,k,INDEX_Energy2) = Rebuilt[INDEX_Energy2];
#else
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        U(i,j,k,INDEX_Energy1+i_Phase) = Rebuilt[INDEX_Energy1+i_Phase];
    }
    U(i,j,k,INDEX_Total_Energy) = Rebuilt[INDEX_Total_Energy];
#endif
    return 0;
}

#endif
#endif
