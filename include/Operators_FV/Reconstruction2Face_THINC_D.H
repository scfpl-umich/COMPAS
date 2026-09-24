// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_THINC_D_H_
#define Reconstruction2Face_THINC_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Reconstruction2Face_THINC_D(amrex::Real beta,
                            amrex::Real sigma,
                            amrex::Real Alpha_bar,
                            int ID_Face)
{
    
    amrex::Real a = (ID_Face == 1) ? 0.0 : 1.0; 
    amrex::Real b = (ID_Face == 1) ? 1.0 : -1.0; 

    amrex::Real beta_ia = beta;

    amrex::Real arg_ln = (std::exp(2.0 * beta_ia * Alpha_bar) - 1.0) / (std::exp(2.0 * beta_ia) - std::exp(2.0 * beta_ia * Alpha_bar));

    amrex::Real c_ia = 0.5/beta_ia * std::log(arg_ln);

    amrex::Real arg_tanh = beta_ia * (a + 0.5 * b * (sigma + 1.0) + c_ia);

    amrex::Real Alpha_face = 0.5 * (1.0 + std::tanh(arg_tanh));

    Alpha_face = amrex::min(1.0, amrex::max(0.0, Alpha_face));

    return Alpha_face;
}

#endif