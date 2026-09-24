// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Mapping_Backward_H_
#define Mapping_Backward_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Mapping_Backward(amrex::Real& Psi,
                      amrex::Real const& PF,
                      amrex::Real const& Eta,
                      int ID_Derivative = 0)
{
    amrex::Real ep = 1.0e-12;
    amrex::Real PF_temp = PF;
    if (PF_temp < ep){
        PF_temp = ep;
    }
    else if (PF_temp > 1.0 - ep){
        PF_temp = 1.0 - ep;
    }
    
    if (ID_Derivative == 0){
        Psi = std::sqrt(2.0)*Eta*std::atanh(2.0*PF_temp - 1.0);
    }
}

#endif