// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseFieldMechanism_CAC_H_
#define PhaseFieldMechanism_CAC_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseFieldMechanism_CAC(amrex::Real& PF_Mechanism,
                             amrex::Real const& PF_Mechanism_AC,
                             amrex::Real const& Integral_Mechanism_AC,
                             amrex::Real const& PF_Weight,
                             amrex::Real const& Integral_Weight)
{
    amrex::Real LagrangeMultiplier;
    if (Integral_Weight == 0.0){
        LagrangeMultiplier = 0.0;
    }
    else{
        LagrangeMultiplier = -Integral_Mechanism_AC/Integral_Weight;
    }
    PF_Mechanism = PF_Mechanism_AC + PF_Weight*LagrangeMultiplier;
}

#endif