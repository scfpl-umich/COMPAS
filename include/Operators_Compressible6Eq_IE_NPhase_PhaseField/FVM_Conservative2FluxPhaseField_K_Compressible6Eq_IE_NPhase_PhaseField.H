// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxPhaseField_K_Compressible6Eq_IE_NPhase_PhaseField_H_
#define FVM_Conservative2FluxPhaseField_K_Compressible6Eq_IE_NPhase_PhaseField_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE

void FVM_Conservative2FluxPhaseField_K(int i, int j, int k,
                                       amrex::Array<amrex::Real, NSTATE>& Flux,
                                       amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
                                       amrex::Real& cmax,
                                       amrex::Array4<amrex::Real const> const& Conservative,
                                       amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                                       Parm const& EOS_Parameter,
                                       amrex::Real const& dx,
                                       int ID_Dimension,
                                       FiniteVolume_Parameter const& FiniteVolume_Parm,
                                       PhaseField_Parameter const& PhaseField_Parm)
{
    // int ID_Bound   = FiniteVolume_Parm.ID_Bound;
    int ID_WENO    = FiniteVolume_Parm.ID_WENO;
    int ID_Limiter = FiniteVolume_Parm.ID_Limiter;
    int ID_Poly    = FiniteVolume_Parm.ID_Poly;

    // Reconstruction
    amrex::Array<amrex::Real, NSTATE> ConservativeL;
    amrex::Array<amrex::Real, NSTATE> ConservativeR;

    FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                             ConservativeL, ConservativeR,
                                             Conservative,
                                             EOS_Parameter,
                                             dx,
                                             ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);
    
    //Flux_PhaseField
    int ID_Flux           = PhaseField_Parm.ID_Flux;
    int ID_FluxPhase2Flux = PhaseField_Parm.ID_FluxPhase2Flux;

    int ID_BE             = 2;  // PhaseField_Parm.ID_Formulation;
    
    FluxPhaseField(Flux,
                   cmax,
                   ConservativeL,
                   ConservativeR,
                   Grad_Q,
                   EOS_Parameter,
                   ID_Flux,
                   ID_FluxPhase2Flux,
                   ID_BE
                );

    // The Phase-Field step has no non-conservative face terms in this model
    for (int iState = 0; iState < NC_TERMS; iState++){
        FluxNC[iState] = 0.0;
    }
    
}

#endif
#endif
