// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Primitive_Compressible6Eq_H_
#define Conservative2Primitive_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Primitive(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Primitive,
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
    amrex::Real Alpha2;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       Alpha1_rho1_E1,
                       Alpha2_rho2_E2,
                       Alpha1,
                       Conservative);
    Alpha2 = 1.0 - Alpha1;
    
    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;
    AMREX_D_TERM( amrex::Real u = rho_u/rho;,
                  amrex::Real v = rho_v/rho;,
                  amrex::Real w = rho_w/rho );
    amrex::Real Alpha1_rho1_e1 = Alpha1_rho1_E1 - (0.5)*Alpha1_rho1*( AMREX_D_TERM( u*u,
                                                                                   +v*v,
                                                                                   +w*w ) );
    amrex::Real Alpha2_rho2_e2 = Alpha2_rho2_E2 - (0.5)*Alpha2_rho2*( AMREX_D_TERM( u*u,
                                                                                   +v*v,
                                                                                   +w*w ) );
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = false;//From Alpha_rho_e to Alpha_p
    int ID_Phase;
    ID_Phase = 1;
    amrex::Real Alpha1_p1;
    EquationOfState_Phase(Alpha1_p1,
                          Alpha1_rho1_e1,
                          Alpha1,
                          Alpha1_rho1,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    ID_Phase = 2;
    amrex::Real Alpha2_p2;
    EquationOfState_Phase(Alpha2_p2,
                          Alpha2_rho2_e2,
                          Alpha2,
                          Alpha2_rho2,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    
    State2Primitive(Primitive,
                    Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    Alpha1_p1,
                    Alpha2_p2,
                    Alpha1);
    
}

#endif
#endif