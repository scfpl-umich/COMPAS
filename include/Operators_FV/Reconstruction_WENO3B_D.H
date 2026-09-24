// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_WENO3B_D_H_
#define Reconstruction_WENO3B_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_WENO3B_D(amrex::Real& f_Reconstruction,
                             amrex::Real const& x,
                             amrex::Real const& f_m,
                             amrex::Real const& f_c,
                             amrex::Real const& f_p,
                             amrex::Real const& f_max,
                             amrex::Real const& f_min,
                             amrex::Real const& x_c,
                             amrex::Real const& dx,
                             int ID_BoundType)
{
    int ID_Face;
    
    //fL_fp
    ID_Face = 1;
    amrex::Real fL_fp;
    Reconstruction2Face_WENO3B_D(fL_fp,
                                 f_m,
                                 f_c,
                                 f_p,
                                 f_max,
                                 f_min,
                                 ID_Face,
                                 ID_BoundType);
    
    //fR_fm
    ID_Face = -1;
    amrex::Real fR_fm;
    Reconstruction2Face_WENO3B_D(fR_fm,
                                 f_m,
                                 f_c,
                                 f_p,
                                 f_max,
                                 f_min,
                                 ID_Face,
                                 ID_BoundType);
    
    Reconstruction_CenterFace2(f_Reconstruction,
                               x,
                               f_c,
                               x_c,
                               fL_fp,
                               fR_fm,
                               dx);
    
}

#endif