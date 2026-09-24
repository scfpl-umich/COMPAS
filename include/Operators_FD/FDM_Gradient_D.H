// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Gradient_D_H_
#define FDM_Gradient_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_D(amrex::Real& Grad_f,
                    amrex::Real const& f_fm,
                    amrex::Real const& f_fp,
                    amrex::Real const& dx)
{
    Grad_f = (f_fp - f_fm)/dx;
}

#endif