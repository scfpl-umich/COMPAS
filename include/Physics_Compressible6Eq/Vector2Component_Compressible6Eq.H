// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Vector2Component_Compressible6Eq_H_
#define Vector2Component_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Vector2Component(amrex::Real& Mass1,
                      amrex::Real& Mass2,
         AMREX_D_DECL(amrex::Real& Momentum_x,
                      amrex::Real& Momentum_y,
                      amrex::Real& Momentum_z),
                      amrex::Real& Energy1,
                      amrex::Real& Energy2,
                      amrex::Real& VolumeFraction1,
                      amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Vector)
{
    Mass1           = Vector[0];
    Mass2           = Vector[1];
    Momentum_x      = Vector[2];
#if (AMREX_SPACEDIM > 1)
    Momentum_y      = Vector[3];
#endif
#if (AMREX_SPACEDIM > 2)
    Momentum_z      = Vector[4];
#endif
    Energy1         = Vector[2+AMREX_SPACEDIM];
    Energy2         = Vector[3+AMREX_SPACEDIM];
    VolumeFraction1 = Vector[4+AMREX_SPACEDIM];
}

#endif
#endif