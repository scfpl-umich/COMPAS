// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativePhase2Conservative_Compressible6Eq_H_
#define ConservativePhase2Conservative_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativePhase2Conservative(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Conservative,
                                    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhase)
{
    VectorPhase2Vector(Conservative,
                       ConservativePhase);
    
}

#endif
#endif