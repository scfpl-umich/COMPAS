// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VectorPhase2Component_Compressible_TwoPhase_H_
#define VectorPhase2Component_Compressible_TwoPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VectorPhase2Component(amrex::Real& Mass1,
                           amrex::Real& Mass2,
                           amrex::Real& Momentum1_x,
                           amrex::Real& Momentum2_x,
#if (AMREX_SPACEDIM > 1)
                           amrex::Real& Momentum1_y,
                           amrex::Real& Momentum2_y,
#endif
#if (AMREX_SPACEDIM > 2)
                           amrex::Real& Momentum1_z,
                           amrex::Real& Momentum2_z,
#endif
                           amrex::Real& Energy1,
                           amrex::Real& Energy2,
                           amrex::Real& VolumeFraction1,
                           amrex::Real& VolumeFraction2,
                           amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& VectorPhase)
{
    Mass1           = VectorPhase[0];
    Mass2           = VectorPhase[1];
    Momentum1_x     = VectorPhase[2];
    Momentum2_x     = VectorPhase[3];
#if (AMREX_SPACEDIM > 1)
    Momentum1_y     = VectorPhase[4];
    Momentum2_y     = VectorPhase[5];
#endif
#if (AMREX_SPACEDIM > 2)
    Momentum1_z     = VectorPhase[6];
    Momentum2_z     = VectorPhase[7];
#endif
    Energy1         = VectorPhase[(1+AMREX_SPACEDIM)*2];
    Energy2         = VectorPhase[(1+AMREX_SPACEDIM)*2+1];
    VolumeFraction1 = VectorPhase[(2+AMREX_SPACEDIM)*2];
    VolumeFraction2 = VectorPhase[(2+AMREX_SPACEDIM)*2+1];
}

#endif
#endif