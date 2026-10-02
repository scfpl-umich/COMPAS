// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Viscosity_Compressible6Eq_H_
#define Conservative2Viscosity_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

// Mixture viscosities of a cell, mu = alpha_1 mu_1 + alpha_2 mu_2 and mub = alpha_1 muB_1 + alpha_2 muB_2
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Viscosity(amrex::Real& mu,
                            amrex::Real& mub,
                            amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                            Parm const& EOS_Parameter)
{

    amrex::Real A1 = Conservative[INDEX_VolumeFraction1];
    amrex::Real A2 = 1. - A1;

    mu = A1*EOS_Parameter.mu_1 + A2*EOS_Parameter.mu_2;
    mub = A1*EOS_Parameter.muB_1 + A2*EOS_Parameter.muB_2;
}

#endif
#endif
