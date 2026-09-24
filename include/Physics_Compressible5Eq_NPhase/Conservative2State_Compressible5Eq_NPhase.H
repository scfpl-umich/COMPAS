// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2State_Compressible5Eq_NPhase_H_
#define Conservative2State_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2State(amrex::Array<amrex::Real,NPHASE>& Alpha_rho,
           AMREX_D_DECL(amrex::Real& rho_u,
                        amrex::Real& rho_v,
                        amrex::Real& rho_w),
                        amrex::Real& rho_E,
                        amrex::Array<amrex::Real,NPHASE>& Alpha,
                        amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative)
{
    Vector2Component(Alpha_rho,
        AMREX_D_DECL(rho_u, rho_v, rho_w),
                     rho_E,
                     Alpha,
                     Conservative);
}

#endif
#endif