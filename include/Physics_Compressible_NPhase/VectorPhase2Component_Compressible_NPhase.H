// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VectorPhase2Component_Compressible_NPhase_H_
#define VectorPhase2Component_Compressible_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VectorPhase2Component(amrex::Array<amrex::Real,NPHASE>& Mass,
                           amrex::Array<amrex::Real,NPHASE>& Momentum_x,
#if (AMREX_SPACEDIM > 1)
                           amrex::Array<amrex::Real,NPHASE>& Momentum_y,
#endif
#if (AMREX_SPACEDIM > 2)
                           amrex::Array<amrex::Real,NPHASE>& Momentum_z,
#endif
                           amrex::Array<amrex::Real,NPHASE>& Energy,
                           amrex::Array<amrex::Real,NPHASE>& VolumeFraction,
                           amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& VectorPhase)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Mass          [i_Phase] = VectorPhase[i_Phase];
        Momentum_x    [i_Phase] = VectorPhase[NPHASE+i_Phase];
#if (AMREX_SPACEDIM > 1)
        Momentum_y    [i_Phase] = VectorPhase[2*NPHASE+i_Phase];
#endif
#if (AMREX_SPACEDIM > 2)
        Momentum_z    [i_Phase] = VectorPhase[3*NPHASE+i_Phase];
#endif
        Energy        [i_Phase] = VectorPhase[(1+AMREX_SPACEDIM)*NPHASE+i_Phase];
        VolumeFraction[i_Phase] = VectorPhase[(2+AMREX_SPACEDIM)*NPHASE+i_Phase];
    }
}

#endif
#endif