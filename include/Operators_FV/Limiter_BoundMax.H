// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_BoundMax_H_
#define Limiter_BoundMax_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_BoundMax(amrex::Real& Theta,
                      amrex::Real const& f_max,
                      amrex::Real const& f_c,
                      amrex::Real const& fC,
                      amrex::Real const& fL_fp,
                      amrex::Real const& fR_fm)
{
    amrex::Real max_f = std::max({ fC,fL_fp,fR_fm });
    Theta = std::max(0.0, f_max - f_c)/std::max(0.0, max_f - f_c);
    Theta = std::min( 1.0,Theta );
}

#endif