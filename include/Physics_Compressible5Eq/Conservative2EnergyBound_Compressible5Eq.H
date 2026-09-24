// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2EnergyBound_Compressible5Eq_H_
#define Conservative2EnergyBound_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2EnergyBound(amrex::Real& EnergyBound,
                              amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                              Parm const& EOS_Parameter)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(amrex::Real rho_u;,
                 amrex::Real rho_v;,
                 amrex::Real rho_w );
    amrex::Real rho_E;
    amrex::Real Alpha1;
    amrex::Real Alpha2;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha1,
                       Conservative);
    Alpha2 = 1.0 - Alpha1;
    
    amrex::Real rho = Alpha1_rho1
                    + Alpha2_rho2;
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Pressure to internal energy
    int ID_Phase;
    
    ID_Phase = 1;//Phase 1
    amrex::Real p_min1 = EOS_Parameter.p_min[0];
    amrex::Real Alpha1_rho1_e1_min;
    EquationOfState_Phase(Alpha1_rho1_e1_min,
                          p_min1,
                          Alpha1,
                          Alpha1_rho1,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    ID_Phase = 2;//Phase 2
    amrex::Real p_min2 = EOS_Parameter.p_min[1];
    amrex::Real Alpha2_rho2_e2_min;
    EquationOfState_Phase(Alpha2_rho2_e2_min,
                          p_min2,
                          Alpha2,
                          Alpha2_rho2,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    amrex::Real rho_e_min = Alpha1_rho1_e1_min
                          + Alpha2_rho2_e2_min;
    EnergyBound = rho_E - (0.5)*( AMREX_D_TERM(rho_u*rho_u,+rho_v*rho_v,+rho_w*rho_w) )/rho - rho_e_min;
    
}

#endif
#endif