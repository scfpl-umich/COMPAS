// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef rho_c2_Phase_Compressible6Eq_H_
#define rho_c2_Phase_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void rho_c2_Phase(amrex::Real& Alpha_rho_c2,
                  amrex::Real const& Alpha_p,
                  amrex::Real const& Alpha,
                  amrex::Real const& EOS_AP,
                  amrex::Real const& EOS_BP)
{
    Alpha_rho_c2 = ( (EOS_AP + 1.0)*Alpha_p + Alpha*EOS_BP )/EOS_AP;
}

#endif
#endif