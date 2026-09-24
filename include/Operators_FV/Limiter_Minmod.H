// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_Minmod_H_
#define Limiter_Minmod_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_Minmod(amrex::Real& Output,
                    amrex::Real const& Input_a,
                    amrex::Real const& Input_b)
{
    Output = std::min( std::max(Input_a,0.0),std::max(Input_b,0.0) )
            +std::max( std::min(Input_a,0.0),std::min(Input_b,0.0) );
}
#endif