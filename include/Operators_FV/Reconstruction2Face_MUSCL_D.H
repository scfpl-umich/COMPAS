// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_MUSCL_D_H_
#define Reconstruction2Face_MUSCL_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_MUSCL_D(amrex::Real& f_Face,
                                 amrex::Real const& f_m,
                                 amrex::Real const& f_c,
                                 amrex::Real const& f_p,
                                 amrex::Real const& dx,
                                 int ID_Face,
                                 int ID_Limiter)
{
    amrex::Real Slope;
    amrex::Real Slope_m = (f_c - f_m) / dx;
    amrex::Real Slope_c = (f_p - f_m) / (2.0*dx);
    amrex::Real Slope_p = (f_p - f_c) / dx;
    Slope_TVD(Slope,
              Slope_m,
              Slope_c,
              Slope_p,
              ID_Limiter);

    Reconstruction2Face_Linear(f_Face,
                               f_c,
                               Slope,
                               dx,
                               ID_Face);
    
}

#endif