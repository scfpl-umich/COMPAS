// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_WENO3_D_H_
#define Reconstruction2Face_WENO3_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_WENO3_D(amrex::Real& f_Face,
                                 amrex::Real const& f_m,
                                 amrex::Real const& f_c,
                                 amrex::Real const& f_p,
                                 int ID_Face)
{
    const int Num_Stencil = 2;
    amrex::Array<amrex::Real,Num_Stencil> f_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Weight_WENO;
    WENO3(f_WENO,
          Weight_WENO,
          f_m,
          f_c,
          f_p,
          ID_Face);
    
    Reconstruction2Face_WENO(f_Face,
                             f_WENO,
                             Weight_WENO,
                             Num_Stencil);
    
}

#endif