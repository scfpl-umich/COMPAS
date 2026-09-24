// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Velocity_Compressible6Eq_IE_NPhase_H_
#define Conservative2Velocity_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Velocity(AMREX_D_DECL(amrex::Real& u,
                                        amrex::Real& v,
                                        amrex::Real& w),
                           amrex::Array<amrex::Real,NSTATE> const& Conservative)
{
    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;


    AMREX_D_TERM(amrex::Real rho_u;,
                amrex::Real rho_v;,
                amrex::Real rho_w);

    amrex::Real rho = 0.0;
    amrex::Real rho_E;

    Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E, Conservative);

    for(int i_phase=0; i_phase<NPHASE; ++i_phase){
        rho += Alpha_rho[i_phase];
    }
    
    AMREX_D_TERM( u = rho_u/rho;,
                    v = rho_v/rho;,
                    w = rho_w/rho );
}

#endif
#endif