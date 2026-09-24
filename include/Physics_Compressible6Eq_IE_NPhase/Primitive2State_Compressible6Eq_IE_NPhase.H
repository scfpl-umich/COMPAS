// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2State_Compressible6Eq_IE_NPhase_H_
#define Primitive2State_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2State(
        amrex::Array<amrex::Real,NPHASE> & Alpha,
        amrex::Array<amrex::Real,NPHASE> & Alpha_rho,
        amrex::Array<amrex::Real,NPHASE> & Alpha_Ps,
        AMREX_D_DECL(amrex::Real & u,
                     amrex::Real & v,
                     amrex::Real & w),
        amrex::Real & p,
        amrex::Array<amrex::Real,NSTATE> const& Primitive)
{
    p = Primitive[NSTATE-1];

    for(int i_phase=0; i_phase<NPHASE; ++i_phase){
        Alpha[i_phase] = Primitive[i_phase];
        Alpha_rho[i_phase] = Primitive[i_phase + NPHASE];
        Alpha_Ps[i_phase] = Primitive[i_phase + 2*NPHASE];
    }

    u = Primitive[3*NPHASE];
    #if (AMREX_SPACEDIM == 2)
        v = Primitive[3*NPHASE +1];

    #endif
    #if (AMREX_SPACEDIM == 3)
        v = Primitive[3*NPHASE + 1];
        w = Primitive[3*NPHASE + 2];
    #endif
}

#endif
#endif