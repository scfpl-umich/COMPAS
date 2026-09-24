// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VelocityPhaseField_Base_Compressible_NPhase_PhaseField_H_
#define VelocityPhaseField_Base_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VelocityPhaseField_Base(amrex::Array<amrex::Real,NPHASE>& VPhase,
                             amrex::Array<amrex::Real,NPHASE> const& Alpha,
                             amrex::Array<amrex::Real,NPHASE> const& Grad_Q)
{
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real Grad_Qi = Grad_Q[i_Phase];
        
        amrex::Real Sum = 0.0;
        for (int j_Phase = 0; j_Phase < NPHASE; j_Phase++){
            amrex::Real Grad_Qj = Grad_Q[j_Phase];
            amrex::Real Alphaj  = Alpha [j_Phase];
            
            amrex::Real Grad_Qij = Grad_Qi - Grad_Qj;
            Sum += Alphaj*Grad_Qij;
        }
        VPhase[i_Phase] = Sum;
    }
    
}

#endif