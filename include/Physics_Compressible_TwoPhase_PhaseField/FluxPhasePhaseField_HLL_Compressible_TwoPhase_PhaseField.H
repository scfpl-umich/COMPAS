// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_HLL_Compressible_TwoPhase_PhaseField_H_
#define FluxPhasePhaseField_HLL_Compressible_TwoPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_HLL(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& FluxPhase,
                             amrex::Real& cmax,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseL,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhaseR,
                             amrex::Real const& Grad_Q,
                             Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> FluxPhaseL;
    amrex::Real VPhase1L;
    amrex::Real VPhase2L;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseL,
                                          VPhase1L, VPhase2L,
                                          ConservativePhaseL,
                                          Grad_Q,
                                          EOS_Parameter);
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> FluxPhaseR;
    amrex::Real VPhase1R;
    amrex::Real VPhase2R;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseR,
                                          VPhase1R, VPhase2R,
                                          ConservativePhaseR,
                                          Grad_Q,
                                          EOS_Parameter);
    
    amrex::Real Sp = std::max({ 0.0, -VPhase1L, -VPhase2L, -VPhase1R, -VPhase2R }) + 1.0e-10;
    amrex::Real Sm = std::min({ 0.0, -VPhase1L, -VPhase2L, -VPhase1R, -VPhase2R }) - 1.0e-10;
    
    for (int i_State = 0; i_State < (3+AMREX_SPACEDIM)*2; i_State++){
        FluxPhase[i_State] = ( Sp*FluxPhaseL[i_State] - Sm*FluxPhaseR[i_State] - Sp*Sm*(ConservativePhaseR[i_State] - ConservativePhaseL[i_State]) )/(Sp - Sm);
    }
    
    cmax = std::max( -Sm, Sp );
    
}

#endif