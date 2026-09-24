// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseFieldMechanism_AC_H_
#define PhaseFieldMechanism_AC_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseFieldMechanism_AC(amrex::Real& PF_Mechanism,
                            amrex::Real const& ChemicalPotential,
                            amrex::Real const& Mob)
{
    PF_Mechanism = -Mob*ChemicalPotential;
}

#endif