// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2EnergyBound_Compressible6Eq_IE_NPhase_H_
#define Conservative2EnergyBound_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2EnergyBound(amrex::Real& EnergyBound,
                              amrex::Array<amrex::Real,NSTATE> const& Conservative,
                              Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e;
    AMREX_D_TERM(amrex::Real rho_u;,
                 amrex::Real rho_v;,
                 amrex::Real rho_w );
    amrex::Real rho_E;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Conservative2State(Alpha,
                       Alpha_rho,
                       Alpha_rho_e,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Conservative);
    
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Pressure to Internal Energy
    amrex::Real p_mini;
    amrex::Real Alphai_rhoi_ei_min;
    amrex::Real rho_e_min = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        p_mini = EOS_Parameter.p_min[i_Phase];
        EquationOfState_Phase(Alphai_rhoi_ei_min,
                              p_mini,
                              Alpha[i_Phase],
                              Alpha_rho[i_Phase],
                              EOS_Parameter,
                              i_Phase,
                              PressureBase,
                              OutputInternalEnergy);
        rho_e_min += Alphai_rhoi_ei_min;
    }
    
    EnergyBound = rho_E - (0.5)*( AMREX_D_TERM(rho_u*rho_u,+rho_v*rho_v,+rho_w*rho_w) )/rho - rho_e_min;
    
}

#endif
#endif