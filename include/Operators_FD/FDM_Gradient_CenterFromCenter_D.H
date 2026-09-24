// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Gradient_CenterFromCenter_D_H_
#define FDM_Gradient_CenterFromCenter_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Gradient_CenterFromCenter_D(amrex::Real& Grad_f,
                                     amrex::Real const& f_m,
                                     amrex::Real const& f_c,
                                     amrex::Real const& f_p,
                                     amrex::Real const& dx)
{
    amrex::Real Grad_fm;
    FDM_Gradient_D(Grad_fm,
                   f_m,
                   f_c,
                   dx);
    
    amrex::Real Grad_fp;
    FDM_Gradient_D(Grad_fp,
                   f_c,
                   f_p,
                   dx);
    
    Interpolation2Face_Linear_D(Grad_f,
                                Grad_fm,
                                Grad_fp);
    
}

#endif