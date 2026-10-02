// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2PhaseTemperature_Compressible6Eq_H_
#define Conservative2PhaseTemperature_Compressible6Eq_H_

#include<Macros.H>
#include<PhaseDiffusion_D.H>
#if (PHYSICS==SIXEQS)

// =======================================
// Temperature of each phase from its own EOS,
//   T_k = (alpha_k rho_k e_k - alpha_k rho_k D_k - alpha_k B_k/(1 + A_k))/(alpha_k rho_k cv_k),
// with alpha_k rho_k e_k = alpha_k rho_k E_k - alpha_k rho_k |u|^2/2. Present[k] tells where
// phase k is present (PhaseDiffusion_Present); T_k is set to 0 where it is not.
// Index 0 is phase 1 and index 1 is phase 2.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2PhaseTemperature(amrex::Array<amrex::Real,2>& Temperature,
                                   amrex::Array<bool,2>& Present,
                                   amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                   Parm const& EOS_Parameter)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM( amrex::Real rho_u;,
                  amrex::Real rho_v;,
                  amrex::Real rho_w );
    amrex::Real Alpha1_rho1_E1;
    amrex::Real Alpha2_rho2_E2;
    amrex::Real Alpha1;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       Alpha1_rho1_E1,
                       Alpha2_rho2_E2,
                       Alpha1,
                       Conservative);
    amrex::Real Alpha2 = 1.0 - Alpha1;

    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;
    AMREX_D_TERM( amrex::Real u = rho_u/rho;,
                  amrex::Real v = rho_v/rho;,
                  amrex::Real w = rho_w/rho );
    amrex::Real Vel2 = AMREX_D_TERM( u*u, +v*v, +w*w );

    amrex::Array<amrex::Real,2> const Alpha     = {Alpha1, Alpha2};
    amrex::Array<amrex::Real,2> const Alpha_rho = {Alpha1_rho1, Alpha2_rho2};
    amrex::Array<amrex::Real,2> const Alpha_rho_E = {Alpha1_rho1_E1, Alpha2_rho2_E2};

    bool PressureBase = false;//Temperature
    bool OutputInternalEnergy = false;//From Alpha_rho_e to Alpha_rho_T
    for (int k = 0; k < 2; ++k){
        Present[k] = PhaseDiffusion_Present(Alpha[k], Alpha_rho[k]);
        Temperature[k] = 0.0;
        if (Present[k]){
            amrex::Real Alpha_rho_e = Alpha_rho_E[k] - (0.5)*Alpha_rho[k]*Vel2;
            amrex::Real Alpha_rho_T;
            EquationOfState_Phase(Alpha_rho_T,
                                  Alpha_rho_e,
                                  Alpha[k],
                                  Alpha_rho[k],
                                  EOS_Parameter,
                                  k+1,
                                  PressureBase,
                                  OutputInternalEnergy);
            Temperature[k] = Alpha_rho_T/Alpha_rho[k];
        }
    }
}

#endif
#endif
