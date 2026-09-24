// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef State2Conservative_Compressible6Eq_H_
#define State2Conservative_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void State2Conservative(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Conservative,
                        amrex::Real const& Alpha1_rho1,
                        amrex::Real const& Alpha2_rho2,
           AMREX_D_DECL(amrex::Real const& rho_u,
                        amrex::Real const& rho_v,
                        amrex::Real const& rho_w),
                        amrex::Real const& Alpha1_rho1_E1,
                        amrex::Real const& Alpha2_rho2_E2,
                        amrex::Real const& Alpha1)
{
    Component2Vector(Conservative,
                     Alpha1_rho1,
                     Alpha2_rho2,
        AMREX_D_DECL(rho_u, rho_v, rho_w),
                     Alpha1_rho1_E1,
                     Alpha2_rho2_E2,
                     Alpha1);
}

#endif
#endif