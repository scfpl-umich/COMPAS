// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Vertex_Linear_D_H_
#define Interpolation2Vertex_Linear_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Vertex_Linear_D(amrex::Real& f_Vertex,
                                   amrex::Real const& f_c,
                                   amrex::Real const& f_1m,
                                   amrex::Real const& f_2m,
                                   amrex::Real const& f_3m,
                                   amrex::Real const& f_1m_2m,
                                   amrex::Real const& f_1m_3m,
                                   amrex::Real const& f_2m_3m,
                                   amrex::Real const& f_1m_2m_3m)
{
    f_Vertex = 0.125*f_c
             + 0.125*f_1m
             + 0.125*f_2m
             + 0.125*f_3m
             + 0.125*f_1m_2m
             + 0.125*f_1m_3m
             + 0.125*f_2m_3m
             + 0.125*f_1m_2m_3m;
}

#endif