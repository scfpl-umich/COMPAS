// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Component2Vector_Compressible6Eq_H_
#define Component2Vector_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Component2Vector(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Vector,
                      amrex::Real const& Mass1,
                      amrex::Real const& Mass2,
         AMREX_D_DECL(amrex::Real const& Momentum_x,
                      amrex::Real const& Momentum_y,
                      amrex::Real const& Momentum_z),
                      amrex::Real const& Energy1,
                      amrex::Real const& Energy2,
                      amrex::Real const& VolumeFraction1)
{
    Vector[0] = Mass1;
    Vector[1] = Mass2;
    Vector[2] = Momentum_x;
#if (AMREX_SPACEDIM > 1)
    Vector[3] = Momentum_y;
#endif
#if (AMREX_SPACEDIM > 2)
    Vector[4] = Momentum_z;
#endif
    Vector[2+AMREX_SPACEDIM] = Energy1;
    Vector[3+AMREX_SPACEDIM] = Energy2;
    Vector[4+AMREX_SPACEDIM] = VolumeFraction1;
}

#endif
#endif