// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Mapping_Forward_H_
#define Mapping_Forward_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Mapping_Forward(amrex::Real& PF,
                     amrex::Real const& Psi,
                     amrex::Real const& Eta,
                     int ID_Derivative = 0)
{
    //PF = F(Psi)
    if (ID_Derivative == 0){
        PF = 0.5 + 0.5*std::tanh( Psi/(std::sqrt(2.0)*Eta) );
    }
    //PF = dF/dPsi
    else if (ID_Derivative == 1){
        PF = 0.5/(std::sqrt(2.0)*Eta)/std::cosh( Psi/(std::sqrt(2.0)*Eta) )/std::cosh( Psi/(std::sqrt(2.0)*Eta) );
    }
    
}

#endif