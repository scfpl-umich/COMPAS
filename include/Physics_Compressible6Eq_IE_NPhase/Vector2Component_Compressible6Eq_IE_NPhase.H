// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Vector2Component_Compressible6Eq_IE_NPhase_H_
#define Vector2Component_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Vector2Component(amrex::Array<amrex::Real,NPHASE>& VolumeFraction,
                      amrex::Array<amrex::Real,NPHASE>& Mass,
                      amrex::Array<amrex::Real,NPHASE>& Energy_Phase,
         AMREX_D_DECL(amrex::Real& Momentum_x,
                      amrex::Real& Momentum_y,
                      amrex::Real& Momentum_z),
                      amrex::Real& Energy,
                      amrex::Array<amrex::Real, NSTATE> const& Vector)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        
        VolumeFraction[i_Phase] = Vector[i_Phase];
        Mass[i_Phase]           = Vector[i_Phase + NPHASE];
        Energy_Phase[i_Phase]   = Vector[i_Phase + 2*NPHASE];

    }
    Momentum_x                  = Vector[3*NPHASE];
#if (AMREX_SPACEDIM > 1)
    Momentum_y                  = Vector[3*NPHASE+1];
#endif
#if (AMREX_SPACEDIM > 2)
    Momentum_z                  = Vector[3*NPHASE+2];
#endif
    Energy                      = Vector[3*NPHASE+AMREX_SPACEDIM];
}

#endif
#endif