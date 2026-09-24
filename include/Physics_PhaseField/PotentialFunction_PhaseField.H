// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PotentialFunction_PhaseField_H_
#define PotentialFunction_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PotentialFunction_PhaseField(amrex::Real& f,
                                  amrex::Real const& x,
                                  int ID_Derivative,
                                  int ID_Potential = 0)
{
    //Double-Well potential
    if (ID_Potential == 0){
        PotentialFunction_DoubleWell(f,
                                     x,
                                     ID_Derivative);
    }
    
}

#endif