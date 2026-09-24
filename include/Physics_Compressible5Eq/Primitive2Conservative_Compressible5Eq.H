// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2Conservative_Compressible5Eq_H_
#define Primitive2Conservative_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2Conservative(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Conservative,
                            amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Primitive,
                            Parm const& EOS_Parameter)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(amrex::Real u;, 
                 amrex::Real v;, 
                 amrex::Real w );
    amrex::Real p;
    amrex::Real Alpha1;
    Primitive2State(Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha1,
                    Primitive);
    
    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;
    AMREX_D_TERM( amrex::Real rho_u = rho*u;,
                  amrex::Real rho_v = rho*v;,
                  amrex::Real rho_w = rho*w );
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Pressure to rho_e
    amrex::Real rho_e;
    EquationOfState_Mixture(rho_e,
                            p,
                            Alpha1,
                            Alpha1_rho1,
                            Alpha2_rho2,
                            EOS_Parameter,
                            PressureBase,
                            OutputInternalEnergy);
    
    amrex::Real rho_E = rho_e + (0.5)*rho*( AMREX_D_TERM( u*u,
                                                         +v*v,
                                                         +w*w ) );
    
    State2Conservative(Conservative,
                       Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha1);
    
}

#endif
#endif