// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef State2Primitive_Compressible5Eq_H_
#define State2Primitive_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void State2Primitive(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Primitive,
                     amrex::Real const& Alpha1_rho1,
                     amrex::Real const& Alpha2_rho2,
        AMREX_D_DECL(amrex::Real const& u,
                     amrex::Real const& v,
                     amrex::Real const& w),
                     amrex::Real const& p,
                     amrex::Real const& Alpha1)
{
    Component2Vector(Primitive,
                     Alpha1_rho1,
                     Alpha2_rho2,
        AMREX_D_DECL(u, v, w),
                     p,
                     Alpha1);
}

#endif
#endif