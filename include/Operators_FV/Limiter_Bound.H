// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_Bound_H_
#define Limiter_Bound_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_Bound(amrex::Real& Theta,
                   amrex::Real const& f_max,
                   amrex::Real const& f_min,
                   amrex::Real const& f_c,
                   amrex::Real const& fL_fp,
                   amrex::Real const& fR_fm,
                   amrex::Real const& Weight,
                   int ID_BoundType)
{
    amrex::Real fC;
    fC_QuadratureRule(fC,
                      f_c,
                      fL_fp,
                      fR_fm,
                      Weight);
    
    //BoundMax
    if (ID_BoundType == 1){
        Limiter_BoundMax(Theta,
                         f_max,
                         f_c,
                         fC,
                         fL_fp,
                         fR_fm);
    }
    //BoundMin
    else if (ID_BoundType == -1){
        Limiter_BoundMin(Theta,
                         f_min,
                         f_c,
                         fC,
                         fL_fp,
                         fR_fm);
    }
    //BoundMaxMin
    else if (ID_BoundType == 0){
        Limiter_BoundMaxMin(Theta,
                            f_max,
                            f_min,
                            f_c,
                            fC,
                            fL_fp,
                            fR_fm);
    }
    
}

#endif