// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2State_Compressible6Eq_IE_NPhase_H_
#define Conservative2State_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2State(amrex::Array<amrex::Real, NPHASE> & Alpha,
                        amrex::Array<amrex::Real, NPHASE> & Alpha_rho,
                        amrex::Array<amrex::Real, NPHASE> & Alpha_rho_e,
           AMREX_D_DECL(amrex::Real& rho_u,
                        amrex::Real& rho_v,
                        amrex::Real& rho_w),
                        amrex::Real& rhoE,
                        amrex::Array<amrex::Real, NSTATE> const& Conservative)
{
    rhoE = Conservative[NSTATE-1];
    rho_u = Conservative[3*NPHASE];

    for(int i_Phase=0; i_Phase < NPHASE; ++i_Phase){
        Alpha[i_Phase] = Conservative[i_Phase];
        Alpha_rho[i_Phase] = Conservative[i_Phase + NPHASE];
        Alpha_rho_e[i_Phase] = Conservative[i_Phase + 2*NPHASE];
    }

    #if (AMREX_SPACEDIM == 2)

        rho_v = Conservative[3*NPHASE + 1];

    #endif
    #if (AMREX_SPACEDIM == 3)
        
        rho_v = Conservative[3*NPHASE + 1];
        rho_w = Conservative[3*NPHASE + 2];
    #endif
}

#endif
#endif