// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Component2VectorPhase_Compressible_NPhase_H_
#define Component2VectorPhase_Compressible_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Component2VectorPhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& VectorPhase,
                           amrex::Array<amrex::Real,NPHASE> const& Mass,
                           amrex::Array<amrex::Real,NPHASE> const& Momentum_x,
#if (AMREX_SPACEDIM > 1)
                           amrex::Array<amrex::Real,NPHASE> const& Momentum_y,
#endif
#if (AMREX_SPACEDIM > 2)
                           amrex::Array<amrex::Real,NPHASE> const& Momentum_z,
#endif
                           amrex::Array<amrex::Real,NPHASE> const& Energy,
                           amrex::Array<amrex::Real,NPHASE> const& VolumeFraction)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        VectorPhase[i_Phase] = Mass[i_Phase];
        VectorPhase[NPHASE+i_Phase] = Momentum_x[i_Phase];
#if (AMREX_SPACEDIM > 1)
        VectorPhase[2*NPHASE+i_Phase] = Momentum_y[i_Phase];
#endif
#if (AMREX_SPACEDIM > 2)
        VectorPhase[3*NPHASE+i_Phase] = Momentum_z[i_Phase];
#endif
        VectorPhase[(1+AMREX_SPACEDIM)*NPHASE+i_Phase] = Energy[i_Phase];
        VectorPhase[(2+AMREX_SPACEDIM)*NPHASE+i_Phase] = VolumeFraction[i_Phase];
    }
}

#endif
#endif