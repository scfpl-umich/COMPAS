// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseFieldFlux_CDI_H_
#define PhaseFieldFlux_CDI_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseFieldFlux_CDI(amrex::Real& PF_Flux,
                        amrex::Real const& Grad_PF,
                        amrex::Real const& PF,
                        amrex::Real const& Normal,
                        amrex::Real const& Mob,
                        amrex::Real const& Eta)
{
    amrex::Real F1 = (std::sqrt(2.0)/Eta)*PF*(1.0 - PF);
    
    PF_Flux = Mob*(Grad_PF - F1*Normal);
}

#endif