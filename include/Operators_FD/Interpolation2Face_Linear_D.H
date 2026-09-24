// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Face_Linear_D_H_
#define Interpolation2Face_Linear_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Face_Linear_D(amrex::Real& f_Face,
                                 amrex::Real const& f_m,
                                 amrex::Real const& f_c)
{
    f_Face = 0.5*f_m
           + 0.5*f_c;
}

#endif