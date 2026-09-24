// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef State2Primitive_Compressible5Eq_NPhase_H_
#define State2Primitive_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void State2Primitive(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Primitive,
                     amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
        AMREX_D_DECL(amrex::Real const& u,
                     amrex::Real const& v,
                     amrex::Real const& w),
                     amrex::Real const& p,
                     amrex::Array<amrex::Real,NPHASE> const& Alpha)
{
    Component2Vector(Primitive,
                     Alpha_rho,
        AMREX_D_DECL(u, v, w),
                     p,
                     Alpha);
}

#endif
#endif