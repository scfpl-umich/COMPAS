// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Primitive_Compressible5Eq_NPhase_H_
#define Conservative2Primitive_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Primitive(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Primitive,
                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                            Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real rho_u;, 
                  amrex::Real rho_v;, 
                  amrex::Real rho_w );
    amrex::Real rho_E;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Conservative2State(Alpha_rho,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha,
                       Conservative);
    
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    AMREX_D_TERM( amrex::Real u = rho_u/rho;,
                  amrex::Real v = rho_v/rho;,
                  amrex::Real w = rho_w/rho );
    amrex::Real rho_e = rho_E - (0.5)*rho*( AMREX_D_TERM( u*u,
                                                         +v*v,
                                                         +w*w ) );
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = false;//From rho_e to Pressure
    amrex::Real p;
    EquationOfState_Mixture(p,
                            rho_e,
                            Alpha,
                            Alpha_rho,
                            EOS_Parameter,
                            PressureBase,
                            OutputInternalEnergy);
    
    State2Primitive(Primitive,
                    Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha);
    
}

#endif
#endif