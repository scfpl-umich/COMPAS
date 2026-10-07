// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef rho_c2_Mixture_Compressible6Eq_IE_NPhase_H_
#define rho_c2_Mixture_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// Frozen mixture value: rho c^2 = sum_k alpha_k rho_k c_k^2
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void rho_c2_Mixture(amrex::Real& rho_c2,
                    amrex::Array<amrex::Real,NPHASE> const& Alphap_rhop_cp2)
{
    rho_c2 = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho_c2 += Alphap_rhop_cp2[i_Phase];
    }
}

#endif
#endif
