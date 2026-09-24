// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Interpolation2Face_WENO3_D_H_
#define Interpolation2Face_WENO3_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Interpolation2Face_WENO3_D(amrex::Real& f_Face,
                                amrex::Real const& f_mm,
                                amrex::Real const& f_m,
                                amrex::Real const& f_c,
                                amrex::Real const& f_p,
                                amrex::Real const& V)
{
    if (V > 0.0){
        int ID_Face = 1;
        Reconstruction2Face_WENO3_D(f_Face,
                                    f_mm,
                                    f_m,
                                    f_c,
                                    ID_Face);
    }
    else if (V < 0.0){
        int ID_Face = -1;
        Reconstruction2Face_WENO3_D(f_Face,
                                    f_m,
                                    f_c,
                                    f_p,
                                    ID_Face);
    }
    else{
        Interpolation2Face_Linear_D(f_Face,
                                    f_m,
                                    f_c);
    }
    
}

#endif