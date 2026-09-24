// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Primitive_Compressible5Eq_H_
#define Conservative2Primitive_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Primitive(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Primitive,
                            amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                            Parm const& EOS_Parameter)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM( amrex::Real rho_u;, 
                  amrex::Real rho_v;, 
                  amrex::Real rho_w );
    amrex::Real rho_E;
    amrex::Real Alpha1;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha1,
                       Conservative);
    
    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;
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
                            Alpha1,
                            Alpha1_rho1,
                            Alpha2_rho2,
                            EOS_Parameter,
                            PressureBase,
                            OutputInternalEnergy);
    
    State2Primitive(Primitive,
                    Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha1);
    
}

#endif
#endif