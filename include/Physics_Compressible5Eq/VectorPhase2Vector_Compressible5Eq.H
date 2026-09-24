// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VectorPhase2Vector_Compressible5Eq_H_
#define VectorPhase2Vector_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VectorPhase2Vector(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Vector,
                        amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& VectorPhase)
{
    amrex::Real Mass1;
    amrex::Real Mass2;
    amrex::Real Momentum1_x;
    amrex::Real Momentum2_x;
#if (AMREX_SPACEDIM > 1)
    amrex::Real Momentum1_y;
    amrex::Real Momentum2_y;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Momentum1_z;
    amrex::Real Momentum2_z;
#endif
    amrex::Real Energy1;
    amrex::Real Energy2;
    amrex::Real VolumeFraction1;
    amrex::Real VolumeFraction2;
    VectorPhase2Component(Mass1, Mass2,
                          Momentum1_x, Momentum2_x,
#if (AMREX_SPACEDIM > 1)
                          Momentum1_y, Momentum2_y,
#endif
#if (AMREX_SPACEDIM > 2)
                          Momentum1_z, Momentum2_z,
#endif
                          Energy1, Energy2,
                          VolumeFraction1, VolumeFraction2,
                          VectorPhase);

    Component2Vector(Vector,
                     Mass1, Mass2,
        AMREX_D_DECL(Momentum1_x + Momentum2_x,
                     Momentum1_y + Momentum2_y,
                     Momentum1_z + Momentum2_z),
                     Energy1     + Energy2,
                     VolumeFraction1);
    
}

#endif
#endif