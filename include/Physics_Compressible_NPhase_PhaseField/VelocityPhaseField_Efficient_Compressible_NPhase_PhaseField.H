// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VelocityPhaseField_Efficient_Compressible_NPhase_PhaseField_H_
#define VelocityPhaseField_Efficient_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VelocityPhaseField_Efficient(amrex::Array<amrex::Real,NPHASE>& VPhase,
                                  amrex::Array<amrex::Real,NPHASE> const& Alpha,
                                  amrex::Array<amrex::Real,NPHASE> const& Grad_Q)
{
    amrex::Real J = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        J += Alpha[i_Phase]*(1.0 - Alpha[i_Phase])*Grad_Q[i_Phase];
    }
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        VPhase[i_Phase] = (1.0 - Alpha[i_Phase])*Grad_Q[i_Phase] - J;
    }
    
}

#endif