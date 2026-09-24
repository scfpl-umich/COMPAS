// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_Godunov_H_
#define Reconstruction_Godunov_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_Godunov(amrex::Real& f,
                            amrex::Real const& x,
                            amrex::Real const& f_c)
{
    f = f_c;
}

#endif