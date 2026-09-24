// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2State_Compressible5Eq_NPhase_H_
#define Primitive2State_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2State(amrex::Array<amrex::Real,NPHASE>& Alpha_rho,
        AMREX_D_DECL(amrex::Real& u,
                     amrex::Real& v,
                     amrex::Real& w),
                     amrex::Real& p,
                     amrex::Array<amrex::Real,NPHASE>& Alpha,
                     amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Primitive)
{
    Vector2Component(Alpha_rho,
        AMREX_D_DECL(u, v, w),
                     p,
                     Alpha,
                     Primitive);
}

#endif
#endif