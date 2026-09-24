// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_CenterFace2_H_
#define Reconstruction_CenterFace2_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_CenterFace2(amrex::Real& f,
                                amrex::Real const& x,
                                amrex::Real const& f_c,
                                amrex::Real const& x_c,
                                amrex::Real const& fL_fp,
                                amrex::Real const& fR_fm,
                                amrex::Real const& dx)
{
    amrex::Real a0 = (1.0/4.0) * (6.0*f_c - fL_fp - fR_fm );
    amrex::Real a1 = (1.0/dx) * (fL_fp - fR_fm);
    amrex::Real a2 = (3.0/(dx*dx)) * (fL_fp + fR_fm - 2.0*f_c);
    f = a0 + a1*(x - x_c) + a2*(x - x_c)*(x - x_c);
}

#endif