// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_WENO3B_D_H_
#define Reconstruction2Face_WENO3B_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction2Face_WENO3B_D(amrex::Real& f_Face,
                                  amrex::Real const& f_m,
                                  amrex::Real const& f_c,
                                  amrex::Real const& f_p,
                                  amrex::Real const& f_max,
                                  amrex::Real const& f_min,
                                  int ID_Face,
                                  int ID_BoundType)
{
    int ID_FaceRe;

    //fL_fp
    ID_FaceRe = 1;
    amrex::Real fL_fp;
    Reconstruction2Face_WENO3_D(fL_fp,
                                f_m,
                                f_c,
                                f_p,
                                ID_FaceRe);
    
    //fR_fm
    ID_FaceRe = -1;
    amrex::Real fR_fm;
    Reconstruction2Face_WENO3_D(fR_fm,
                                f_m,
                                f_c,
                                f_p,
                                ID_FaceRe);
    
    //Limiter_Bound
    amrex::Real Weight = 1.0/6.0;
    amrex::Real Theta;
    Limiter_Bound(Theta,
                  f_max,
                  f_min,
                  f_c,
                  fL_fp,
                  fR_fm,
                  Weight,
                  ID_BoundType);
    
    //f_Face
    amrex::Real f;
    //fL_fp
    if (ID_Face == 1){
        f = fL_fp;
    }
    //fR_fm
    else if (ID_Face == -1){
        f = fR_fm;
    }
    f_LimiterBound(f_Face,
                   f,
                   f_c,
                   Theta);
    
}

#endif