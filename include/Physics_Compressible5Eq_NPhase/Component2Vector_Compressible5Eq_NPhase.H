// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Component2Vector_Compressible5Eq_NPhase_H_
#define Component2Vector_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Component2Vector(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Vector,
                      amrex::Array<amrex::Real,NPHASE> const& Mass,
         AMREX_D_DECL(amrex::Real const& Momentum_x,
                      amrex::Real const& Momentum_y,
                      amrex::Real const& Momentum_z),
                      amrex::Real const& Energy,
                      amrex::Array<amrex::Real,NPHASE> const& VolumeFraction)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Vector[i_Phase]                         = Mass[i_Phase];
        Vector[NPHASE+AMREX_SPACEDIM+1+i_Phase] = VolumeFraction[i_Phase];
    }
    Vector[NPHASE]                              = Momentum_x;
#if (AMREX_SPACEDIM > 1)
    Vector[NPHASE+1]                            = Momentum_y;
#endif
#if (AMREX_SPACEDIM > 2)
    Vector[NPHASE+2]                            = Momentum_z;
#endif
    Vector[NPHASE+AMREX_SPACEDIM]               = Energy;
}

#endif
#endif