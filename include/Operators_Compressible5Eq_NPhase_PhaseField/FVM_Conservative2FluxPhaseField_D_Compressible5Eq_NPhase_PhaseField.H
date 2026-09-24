// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxPhaseField_D_Compressible5Eq_NPhase_PhaseField_H_
#define FVM_Conservative2FluxPhaseField_D_Compressible5Eq_NPhase_PhaseField_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

template<typename T_Stencil>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxPhaseField_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                                       amrex::Real& cmax,
                                       T_Stencil const& Conservative_Stencil,
                                       amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                                       Parm const& EOS_Parameter,
                                       amrex::Real const& dx,
                                       int ID_WENO,
                                       int ID_Limiter,
                                       int ID_Poly,
                                       int ID_FluxPhaseField,
                                       int ID_FluxPhase2Flux,
                                       int ID_BE)
{
    //Reconstruction
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeL;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeR;
    FVM_Reconstruction2FaceLR_Conservative_D(ConservativeL, ConservativeR,
                                             Conservative_Stencil,
                                             EOS_Parameter,
                                             dx,
                                             ID_WENO, ID_Limiter, ID_Poly);
    
    //Flux_PhaseField
    FluxPhaseField(Flux, cmax,
                   ConservativeL, ConservativeR, Grad_Q,
                   EOS_Parameter,
                   ID_FluxPhaseField, ID_FluxPhase2Flux, ID_BE);
    
}

#endif
#endif