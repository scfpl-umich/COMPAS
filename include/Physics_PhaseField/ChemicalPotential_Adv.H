// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ChemicalPotential_Adv_H_
#define ChemicalPotential_Adv_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ChemicalPotential_Adv(amrex::Real& ChemicalPotential,
                           amrex::Real const& PF,
                           amrex::Real const& Laplace_PF,
                           amrex::Real const& Abs_Grad_PF,
                           amrex::Real const& Curvature,
                           amrex::Real const& Lambda,
                           amrex::Real const& Eta,
                           int ID_Potential = 0)
{
    int ID_Derivative = 1;
    amrex::Real g1_PF;
    PotentialFunction_PhaseField(g1_PF,
                                 PF,
                                 ID_Derivative,
                                 ID_Potential);
    ChemicalPotential = Lambda*( (1.0/Eta/Eta)*g1_PF - Laplace_PF - Abs_Grad_PF*Curvature );
    
}

#endif