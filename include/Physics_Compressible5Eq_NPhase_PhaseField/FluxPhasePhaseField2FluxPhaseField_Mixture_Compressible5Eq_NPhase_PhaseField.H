// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField2FluxPhaseField_Mixture_Compressible5Eq_NPhase_PhaseField_H_
#define FluxPhasePhaseField2FluxPhaseField_Mixture_Compressible5Eq_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField2FluxPhaseField_Mixture(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                                                amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& FluxPhase,
                                                amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& ConservativeL,
                                                amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& ConservativeR,
                                                Parm const& EOS_Parameter)

{
    amrex::Array<amrex::Real,NPHASE> J_rho;
    amrex::Array<amrex::Real,NPHASE> J_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> J_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> J_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> J_rho_E;
    amrex::Array<amrex::Real,NPHASE> J;
    VectorPhase2Component(J_rho,
                          J_rho_u,
#if (AMREX_SPACEDIM > 1)
                          J_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                          J_rho_w,
#endif
                          J_rho_E,
                          J,
                          FluxPhase);

    FluxPhaseField_Mixture(Flux,
                           ConservativeL,
                           ConservativeR,
                           J_rho,
                           J,
                           EOS_Parameter);
    
}

#endif