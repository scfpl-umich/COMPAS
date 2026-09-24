// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef fC_QuadratureRule_H_
#define fC_QuadratureRule_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fC_QuadratureRule(amrex::Real& fC,
                       amrex::Real const& f_c,
                       amrex::Real const& fL_fp,
                       amrex::Real const& fR_fm,
                       amrex::Real const& Weight)
{
    fC = (f_c - Weight*fL_fp - Weight*fR_fm)/(1.0 - 2.0*Weight);
}

#endif