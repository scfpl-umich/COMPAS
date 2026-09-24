// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Edge_Linear_D_H_
#define Interpolation2Edge_Linear_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Edge_Linear_D(amrex::Real& f_Edge,
                                 amrex::Real const& f_c,
                                 amrex::Real const& f_1m,
                                 amrex::Real const& f_2m,
                                 amrex::Real const& f_1m_2m)
{
    f_Edge = 0.25*f_c
           + 0.25*f_1m
           + 0.25*f_2m
           + 0.25*f_1m_2m;
}

#endif