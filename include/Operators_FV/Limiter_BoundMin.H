// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_BoundMin_H_
#define Limiter_BoundMin_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_BoundMin(amrex::Real& Theta,
                      amrex::Real const& f_min,
                      amrex::Real const& f_c,
                      amrex::Real const& fC,
                      amrex::Real const& fL_fp,
                      amrex::Real const& fR_fm)
{
    amrex::Real min_f = std::min({ fC,fL_fp,fR_fm });
    Theta = std::max(0.0, f_c - f_min)/std::max(0.0, f_c - min_f);
    Theta = std::min( 1.0,Theta );
}

#endif