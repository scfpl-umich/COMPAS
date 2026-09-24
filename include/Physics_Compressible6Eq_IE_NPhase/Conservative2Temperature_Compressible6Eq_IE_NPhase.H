// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Temperature_Compressible6Eq_IE_NPhase_H_
#define Conservative2Temperature_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Temperature(amrex::Real& Temperature,
                              amrex::Array<amrex::Real, NSTATE> const& Conservative,
                              Parm const& EOS_Parameter)
{
    // Not sure if we need one Temperature or multiple phasic temperatures
  Temperature = 0.;
  amrex::Array<amrex::Real, NPHASE> Alpha;
  amrex::Array<amrex::Real, NPHASE> Alpha_rho;
  amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;

  AMREX_D_TERM(amrex::Real rho_u;,
               amrex::Real rho_v;,
               amrex::Real rho_w);

  amrex::Real rho = 0.;
  amrex::Real rhoT = 0.;
  amrex::Real rho_E;

  Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E, Conservative);
  
  
  // int ID_PT = 2;//Temperature
  // int ID_Inverse = 1;//From rho_e to Temperature

  bool PressureBase = false;//Pressure
  bool OutputInternalEnergy = false;//From rho_e to Pressure

  for(int i_k=0; i_k < NPHASE; i_k++){
    amrex::Real akrhokTk = 0.;

    EquationOfState_Phase(akrhokTk,
                  Alpha_rho_e[i_k],
                  Alpha[i_k],
                  Alpha_rho[i_k],
                  EOS_Parameter,
                  i_k,
                  PressureBase,
                  OutputInternalEnergy);
    rho += Alpha_rho[i_k];
    rhoT += akrhokTk; 
  }

  Temperature = rhoT/rho;
    
}

#endif
#endif