// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxPhaseField_K_Compressible6Eq_PhaseField_H_
#define FVM_Conservative2FluxPhaseField_K_Compressible6Eq_PhaseField_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxPhaseField_K(int i, int j, int k,
                                       amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Flux,
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
    int ID_WENO    = FiniteVolume_Parm.ID_WENO;
    int ID_Limiter = FiniteVolume_Parm.ID_Limiter;
    int ID_Poly    = FiniteVolume_Parm.ID_Poly;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> ConservativeL;
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> ConservativeR;
    FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                             ConservativeL, ConservativeR,
                                             Conservative,
                                             EOS_Parameter,
                                             dx,
                                             ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);
    
    //Flux_PhaseField
    int ID_Flux           = PhaseField_Parm.ID_Flux;
    int ID_FluxPhase2Flux = PhaseField_Parm.ID_FluxPhase2Flux;
    amrex::Real J1_rho1;
    amrex::Real J2_rho2;
    amrex::Real J1_rho1_EK1;
    amrex::Real J2_rho2_EK2;
    FluxPhaseField(Flux,
                   J1_rho1, J2_rho2,
                   J1_rho1_EK1, J2_rho2_EK2,
                   cmax,
                   ConservativeL,
                   ConservativeR,
                   Grad_Q[0],
                   EOS_Parameter,
                   ID_Flux,
                   ID_FluxPhase2Flux);

    //FluxNC_PhaseField
    for (int iState = 0; iState < NC_TERMS; iState++){
        FluxNC[iState] = 0.0;
    }
    FluxNC[0] = J1_rho1;
    FluxNC[1] = J2_rho2;
    FluxNC[2] = J1_rho1_EK1;
    FluxNC[3] = J2_rho2_EK2;
    
}

#endif
#endif