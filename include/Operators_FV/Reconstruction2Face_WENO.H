// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_WENO_H_
#define Reconstruction2Face_WENO_H_

template<typename T_WENO>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_WENO(amrex::Real& f_Face,
                              T_WENO const& f_WENO,
                              T_WENO const& Weight_WENO,
                              int Num_Stencil)
{
    f_Face = 0.0;
    for (int i_Stencil=0; i_Stencil<Num_Stencil; ++i_Stencil){
        f_Face += Weight_WENO[i_Stencil]*f_WENO[i_Stencil];
    }
}

#endif