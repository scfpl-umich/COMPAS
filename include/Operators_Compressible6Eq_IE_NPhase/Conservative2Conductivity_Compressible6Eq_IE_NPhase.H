// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Conductivity_Compressible6Eq_IE_NPhase_H_
#define Conservative2Conductivity_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Conductivity(amrex::Real& Kappa,
                               amrex::Array<amrex::Real,NSTATE> const& Conservative,
                               Parm const& EOS_Parameter)
{
    Kappa = 0.;

    for(int i_k = 0; i_k < NPHASE; ++i_k){
        Kappa += Conservative[i_k]*EOS_Parameter.kappa[i_k];
    }
    
}

#endif
#endif