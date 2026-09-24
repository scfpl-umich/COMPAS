// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef State2Primitive_Compressible6Eq_IE_NPhase_H_
#define State2Primitive_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
State2Primitive(
        amrex::Array<amrex::Real,NSTATE>& Primitive,
        amrex::Array<amrex::Real,NPHASE> const& Alpha,
        amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
        amrex::Array<amrex::Real,NPHASE> const& Alpha_Ps,
        AMREX_D_DECL(amrex::Real const& u,
                     amrex::Real const& v,
                     amrex::Real const& w),
        amrex::Real const& p)
{
    Primitive[NSTATE-1] = p;

    for(int i_phase=0; i_phase<NPHASE; ++i_phase){
        Primitive[i_phase] = Alpha[i_phase];
        Primitive[i_phase + NPHASE] = Alpha_rho[i_phase];
        Primitive[i_phase + 2*NPHASE] = Alpha_Ps[i_phase];
    }

    Primitive[3*NPHASE] = u;
    #if (AMREX_SPACEDIM == 2)
        Primitive[3*NPHASE+1] = v;

    #endif
    #if (AMREX_SPACEDIM == 3)
        Primitive[3*NPHASE + 1] = v;
        Primitive[3*NPHASE + 2] = w;
    #endif
}

#endif
#endif