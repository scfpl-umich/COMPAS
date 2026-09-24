// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_Godunov_K_H_
#define Reconstruction2Face_Godunov_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_Godunov_K(int i, int j, int k, int n,
                                   amrex::Real& f_Face,
                                   amrex::Array4<amrex::Real const> const& f)
{
    amrex::Real f_c = f(i,j,k,n);
    Reconstruction2Face_Godunov_D(f_Face,
                                  f_c);
}

#endif