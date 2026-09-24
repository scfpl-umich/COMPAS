// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Conductivity_Compressible5Eq_NPhase_H_
#define Conservative2Conductivity_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Conductivity(amrex::Real& Kappa,
                               amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                               Parm const& EOS_Parameter)
{
    Kappa = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Kappa += Conservative[INDEX_VolumeFraction1+i_Phase]*EOS_Parameter.kappa[i_Phase];
    }
    
}

#endif
#endif