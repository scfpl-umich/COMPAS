// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseFieldFlux_ACDI_H_
#define PhaseFieldFlux_ACDI_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseFieldFlux_ACDI(amrex::Real& PF_Flux,
                         amrex::Real const& Grad_PF,
                         amrex::Real const& Psi,
                         amrex::Real const& Normal,
                         amrex::Real const& Mob,
                         amrex::Real const& Eta)
{
    int ID_Derivative = 1;
    amrex::Real F1;
    Mapping_Forward(F1,
                    Psi,
                    Eta,
                    ID_Derivative);
    PF_Flux = Mob*(Grad_PF - F1*Normal);
}

#endif