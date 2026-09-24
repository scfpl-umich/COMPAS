// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxHyperbolic_D_Compressible5Eq_NPhase_H_
#define FVM_Conservative2FluxHyperbolic_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

template<typename T_Stencil>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxHyperbolic_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                                       amrex::Array<amrex::Real,AMREX_SPACEDIM>& Vel,
                                       amrex::Real& cmax,
                                       T_Stencil const& Conservative_Stencil,
                                       Parm const& EOS_Parameter,
                                       amrex::Real const& dx,
                                       int ID_Dimension,
                                       int ID_WENO,
                                       int ID_Limiter,
                                       int ID_Poly,
                                       int ID_K,
                                       int ID_Riemann)
{
    //Reconstruction
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeL;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> ConservativeR;
    FVM_Reconstruction2FaceLR_Conservative_D(ConservativeL, ConservativeR,
                                             Conservative_Stencil,
                                             EOS_Parameter,
                                             dx,
                                             ID_WENO, ID_Limiter, ID_Poly);
    
    //Flux_Hyperbolic
    RiemannSolver(Flux, Vel, cmax,
                  ConservativeL, ConservativeR,
                  EOS_Parameter,
                  ID_Dimension, ID_K, ID_Riemann);
    
}

#endif
#endif