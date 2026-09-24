// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_MC_H_
#define Limiter_MC_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_MC(amrex::Real& Output,
                amrex::Real const& Input_a,
                amrex::Real const& Input_b,
                amrex::Real const& Input_c)
{
    Output = std::min({ std::max(Input_c,0.0),std::max(2.0*Input_a,0.0),std::max(2.0*Input_b,0.0) })
            +std::max({ std::min(Input_c,0.0),std::min(2.0*Input_a,0.0),std::min(2.0*Input_b,0.0) });
}
#endif