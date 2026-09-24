// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField_HLL_Compressible_NPhase_PhaseField_H_
#define FluxPhasePhaseField_HLL_Compressible_NPhase_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField_HLL(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& FluxPhase,
                             amrex::Real& cmax,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhaseL,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& ConservativePhaseR,
                             amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                             Parm const& EOS_Parameter,
                             int ID_Formulation)
{
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> FluxPhaseL;
    amrex::Array<amrex::Real,NPHASE> VPhaseL;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseL,
                                          VPhaseL,
                                          ConservativePhaseL,
                                          Grad_Q,
                                          EOS_Parameter,
                                          ID_Formulation);
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> FluxPhaseR;
    amrex::Array<amrex::Real,NPHASE> VPhaseR;
    ConservativePhase2FluxPhasePhaseField(FluxPhaseR,
                                          VPhaseR,
                                          ConservativePhaseR,
                                          Grad_Q,
                                          EOS_Parameter,
                                          ID_Formulation);

    amrex::Array<amrex::Real,NPHASE> mVPhaseL;
    amrex::Array<amrex::Real,NPHASE> mVPhaseR;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        mVPhaseL[i_Phase] = -VPhaseL[i_Phase];
        mVPhaseR[i_Phase] = -VPhaseR[i_Phase];
    }
    
    amrex::Real max_mVPhaseL = *std::max_element(std::begin(mVPhaseL), std::end(mVPhaseL));
    amrex::Real max_mVPhaseR = *std::max_element(std::begin(mVPhaseR), std::end(mVPhaseR));
    amrex::Real Sp = std::max({ 0.0, max_mVPhaseL, max_mVPhaseR }) + 1.0e-10;
    
    amrex::Real min_mVPhaseL = *std::min_element(std::begin(mVPhaseL), std::end(mVPhaseL));
    amrex::Real min_mVPhaseR = *std::min_element(std::begin(mVPhaseR), std::end(mVPhaseR));
    amrex::Real Sm = std::min({ 0.0, min_mVPhaseL, min_mVPhaseR }) - 1.0e-10;
    
    for (int i_State = 0; i_State < (3+AMREX_SPACEDIM)*NPHASE; i_State++){
        FluxPhase[i_State] = ( Sp*FluxPhaseL[i_State] - Sm*FluxPhaseR[i_State] - Sp*Sm*(ConservativePhaseR[i_State] - ConservativePhaseL[i_State]) )/(Sp - Sm);
    }
    
    cmax = std::max( -Sm, Sp );
    
}

#endif