// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Conductivity_Compressible5Eq_H_
#define Conservative2Conductivity_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Conductivity(amrex::Real& Kappa,
                               amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                               Parm const& EOS_Parameter)
{
    amrex::Real A1 = Conservative[INDEX_VolumeFraction1];
    amrex::Real A2 = 1. - A1;

    Kappa = A1*EOS_Parameter.kappa_1 + A2*EOS_Parameter.kappa_2;
    
}

#endif
#endif