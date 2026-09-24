// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction_CenterFace4_H_
#define Reconstruction_CenterFace4_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Reconstruction_CenterFace4(amrex::Real& f,
                                amrex::Real const& x,
                                amrex::Real const& f_c,
                                amrex::Real const& x_c,
                                amrex::Real const& f_p,
                                amrex::Real const& f_m,
                                amrex::Real const& fL_fp,
                                amrex::Real const& fR_fm,
                                amrex::Real const& dx)
{
    amrex::Real a0 = (1.0/192.0) * ( (f_p + 298.0*f_c + f_m) - 54.0*(fL_fp + fR_fm) );
    amrex::Real a1 = (-1.0/8.0/dx) * ( (f_p - f_m) + 10.0*(fR_fm - fL_fp) );
    amrex::Real a2 = (-1.0/8.0/std::pow(dx,2)) * ( (f_p + 58.0*f_c + f_m) - 30.0*(fL_fp + fR_fm) );
    amrex::Real a3 = (1.0/2.0/std::pow(dx,3)) * ( (f_p - f_m) + 2.0*(fR_fm - fL_fp) );
    amrex::Real a4 = (5.0/12.0/std::pow(dx,4)) * ( (f_p + 10.0*f_c + f_m) - 6.0*(fL_fp + fR_fm) );

    f = a0 + a1*(x - x_c) + a2*std::pow(x - x_c,2) + a3*std::pow(x - x_c,3) + a4*std::pow(x - x_c,4);
}

#endif