// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Vector2Component_Compressible5Eq_NPhase_H_
#define Vector2Component_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Vector2Component(amrex::Array<amrex::Real,NPHASE>& Mass,
         AMREX_D_DECL(amrex::Real& Momentum_x,
                      amrex::Real& Momentum_y,
                      amrex::Real& Momentum_z),
                      amrex::Real& Energy,
                      amrex::Array<amrex::Real,NPHASE>& VolumeFraction,
                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Vector)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Mass[i_Phase]           = Vector[i_Phase];
        VolumeFraction[i_Phase] = Vector[NPHASE+AMREX_SPACEDIM+1+i_Phase];
    }
    Momentum_x                  = Vector[NPHASE];
#if (AMREX_SPACEDIM > 1)
    Momentum_y                  = Vector[NPHASE+1];
#endif
#if (AMREX_SPACEDIM > 2)
    Momentum_z                  = Vector[NPHASE+2];
#endif
    Energy                      = Vector[NPHASE+AMREX_SPACEDIM];
}

#endif
#endif