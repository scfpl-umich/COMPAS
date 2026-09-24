// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef f_Upwind_H_
#define f_Upwind_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void f_Upwind(amrex::Real& f,
              amrex::Real const& fL,
              amrex::Real const& fR,
              amrex::Real const& V)
{
    if (V > 0.0){
        f = fL;
    }
    else if (V < 0.0){
        f = fR;
    }
    else{
        f = 0.5 * (fL + fR);
    }
}

#endif