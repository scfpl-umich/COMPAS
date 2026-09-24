// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef f_LimiterBound_H_
#define f_LimiterBound_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void f_LimiterBound(amrex::Real& f_Limit,
                    amrex::Real const& f,
                    amrex::Real const& f_c,
                    amrex::Real const& Theta)
{
    f_Limit = f_c + Theta*(f - f_c);
}

#endif