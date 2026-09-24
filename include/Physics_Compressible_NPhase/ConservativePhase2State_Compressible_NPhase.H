// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativePhase2State_Compressible_NPhase_H_
#define ConservativePhase2State_Compressible_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativePhase2State(amrex::Array<amrex::Real,NPHASE>& Alpha_rho,
                             amrex::Array<amrex::Real,NPHASE>& Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                             amrex::Array<amrex::Real,NPHASE>& Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                             amrex::Array<amrex::Real,NPHASE>& Alpha_rho_w,
#endif
                             amrex::Array<amrex::Real,NPHASE>& Alpha_rho_E,
                             amrex::Array<amrex::Real,NPHASE>& Alpha,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhase)
{
    VectorPhase2Component(Alpha_rho,
                          Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                          Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                          Alpha_rho_w,
#endif
                          Alpha_rho_E,
                          Alpha,
                          ConservativePhase);
}

#endif
#endif