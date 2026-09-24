// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_Linear_H_
#define Reconstruction2Face_Linear_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_Linear(amrex::Real& f_Face,
                                amrex::Real const& f_c,
                                amrex::Real const& Slope,
                                amrex::Real const& dx,
                                int ID_Face)
{
    f_Face = f_c + amrex::Real(ID_Face)*0.5*Slope*dx;
}

#endif