// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Component2VectorPhase_Compressible_TwoPhase_H_
#define Component2VectorPhase_Compressible_TwoPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Component2VectorPhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& VectorPhase,
                           amrex::Real const& Mass1,
                           amrex::Real const& Mass2,
                           amrex::Real const& Momentum1_x,
                           amrex::Real const& Momentum2_x,
#if (AMREX_SPACEDIM > 1)
                           amrex::Real const& Momentum1_y,
                           amrex::Real const& Momentum2_y,
#endif
#if (AMREX_SPACEDIM > 2)
                           amrex::Real const& Momentum1_z,
                           amrex::Real const& Momentum2_z,
#endif
                           amrex::Real const& Energy1,
                           amrex::Real const& Energy2,
                           amrex::Real const& VolumeFraction1,
                           amrex::Real const& VolumeFraction2)
{
    VectorPhase[0] = Mass1;
    VectorPhase[1] = Mass2;
    VectorPhase[2] = Momentum1_x;
    VectorPhase[3] = Momentum2_x;
#if (AMREX_SPACEDIM > 1)
    VectorPhase[4] = Momentum1_y;
    VectorPhase[5] = Momentum2_y;
#endif
#if (AMREX_SPACEDIM > 2)
    VectorPhase[6] = Momentum1_z;
    VectorPhase[7] = Momentum2_z;
#endif
    VectorPhase[(1+AMREX_SPACEDIM)*2] = Energy1;
    VectorPhase[(1+AMREX_SPACEDIM)*2+1] = Energy2;
    VectorPhase[(2+AMREX_SPACEDIM)*2] = VolumeFraction1;
    VectorPhase[(2+AMREX_SPACEDIM)*2+1] = VolumeFraction2;
}

#endif
#endif