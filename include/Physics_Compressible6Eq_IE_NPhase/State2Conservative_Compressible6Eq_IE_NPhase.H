// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef State2Conservative_Compressible6Eq_IE_NPhase_H_
#define State2Conservative_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
State2Conservative(
        amrex::Array<amrex::Real,NSTATE>& Conservative,
        amrex::Array<amrex::Real,NPHASE> const& Alpha,
        amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
        amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_e,
        AMREX_D_DECL(amrex::Real const& rho_u,
                     amrex::Real const& rho_v,
                     amrex::Real const& rho_w),
        amrex::Real const& rho_E)
{
    Conservative[NSTATE-1] = rho_E;

    for(int i_phase=0; i_phase<NPHASE; ++i_phase){
        Conservative[i_phase] = Alpha[i_phase];
        Conservative[i_phase + NPHASE] = Alpha_rho[i_phase];
        Conservative[i_phase + 2*NPHASE] = Alpha_rho_e[i_phase];
    }

    Conservative[3*NPHASE] = rho_u;
    #if (AMREX_SPACEDIM == 2)
        Conservative[3*NPHASE+1] = rho_v;

    #endif
    #if (AMREX_SPACEDIM == 3)
        Conservative[3*NPHASE + 1] = rho_v;
        Conservative[3*NPHASE + 2] = rho_w;
    #endif
}

#endif
#endif