// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2PhaseTemperature_Compressible6Eq_IE_NPhase_H_
#define Conservative2PhaseTemperature_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#include<PhaseDiffusion_D.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// =======================================
// Temperature of each phase from its own EOS and its internal energy,
//   T_k = (alpha_k rho_k e_k - alpha_k rho_k D_k - phi_k B_k/(1 + A_k))/(alpha_k rho_k cv_k),
// phi_k = alpha_k - b_k alpha_k rho_k (NASG, D_k = q_k: T_k = (p_k + pinf_k)(1/rho_k - b_k)/((gamma_k - 1) cv_k);
// phi_k = alpha_k without a covolume). Present[k] tells where phase k is present
// (PhaseDiffusion_Present); T_k is set to 0 where it is not.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2PhaseTemperature(amrex::Array<amrex::Real,NPHASE>& Temperature,
                                   amrex::Array<bool,NPHASE>& Present,
                                   amrex::Array<amrex::Real,NSTATE> const& Conservative,
                                   Parm const& EOS_Parameter)
{
    bool PressureBase = false;//Temperature
    bool OutputInternalEnergy = false;//From Alpha_rho_e to Alpha_rho_T
    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        amrex::Real const Alpha       = Conservative[INDEX_VolumeFraction1 + i_Phase];
        amrex::Real const Alpha_rho   = Conservative[INDEX_Mass1 + i_Phase];
        amrex::Real const Alpha_rho_e = Conservative[INDEX_Energy1 + i_Phase];
        Present[i_Phase] = PhaseDiffusion_Present(Alpha, Alpha_rho);
        Temperature[i_Phase] = 0.0;
        if (Present[i_Phase]){
            amrex::Real Alpha_rho_T;
            EquationOfState_Phase(Alpha_rho_T,
                                  Alpha_rho_e,
                                  Alpha,
                                  Alpha_rho,
                                  EOS_Parameter,
                                  i_Phase,
                                  PressureBase,
                                  OutputInternalEnergy);
            Temperature[i_Phase] = Alpha_rho_T/Alpha_rho;
        }
    }
}

#endif
#endif
