// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef rho_c2_Mixture_Compressible6Eq_IE_NPhase_H_
#define rho_c2_Mixture_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void rho_c2_Mixture(amrex::Real& rho_c2,
                    amrex::Real const& Alpha1_rho1_c12,
                    amrex::Real const& Alpha2_rho2_c22)
{

    // rho_c2 = Alpha1_rho1_c12 + Alpha2_rho2_c22;

}

#endif
#endif