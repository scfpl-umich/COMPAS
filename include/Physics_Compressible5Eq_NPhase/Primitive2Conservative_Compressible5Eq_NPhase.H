// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2Conservative_Compressible5Eq_NPhase_H_
#define Primitive2Conservative_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2Conservative(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Conservative,
                            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Primitive,
                            Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM(amrex::Real u;, 
                 amrex::Real v;, 
                 amrex::Real w );
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Primitive2State(Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha,
                    Primitive);
    
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    AMREX_D_TERM( amrex::Real rho_u = rho*u;,
                  amrex::Real rho_v = rho*v;,
                  amrex::Real rho_w = rho*w );
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Pressure to rho_e
    amrex::Real rho_e;
    EquationOfState_Mixture(rho_e,
                            p,
                            Alpha,
                            Alpha_rho,
                            EOS_Parameter,
                            PressureBase,
                            OutputInternalEnergy);
    
    amrex::Real rho_E = rho_e + (0.5)*rho*( AMREX_D_TERM( u*u,
                                                         +v*v,
                                                         +w*w ) );
    
    State2Conservative(Conservative,
                       Alpha_rho,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha);
    
}

#endif
#endif