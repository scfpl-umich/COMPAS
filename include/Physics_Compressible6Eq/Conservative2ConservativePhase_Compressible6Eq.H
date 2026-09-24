// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2ConservativePhase_Compressible6Eq_H_
#define Conservative2ConservativePhase_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2ConservativePhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& ConservativePhase,
                                    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(amrex::Real u;, 
                 amrex::Real v;, 
                 amrex::Real w );
    amrex::Real Alpha1_p1;
    amrex::Real Alpha2_p2;
    amrex::Real Alpha1;
    amrex::Real Alpha2;
    Primitive2State(Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    Alpha1_p1,
                    Alpha2_p2,
                    Alpha1,
                    Primitive);
    Alpha2 = 1.0 - Alpha1;
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Alpha_p to Alpha_rho_e
    int ID_Phase;
    ID_Phase = 1;
    amrex::Real Alpha1_rho1_e1;
    EquationOfState_Phase(Alpha1_rho1_e1,
                          Alpha1_p1,
                          Alpha1,
                          Alpha1_rho1,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    ID_Phase = 2;
    amrex::Real Alpha2_rho2_e2;
    EquationOfState_Phase(Alpha2_rho2_e2,
                          Alpha2_p2,
                          Alpha2,
                          Alpha2_rho2,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
    
    State2ConservativePhase(ConservativePhase,
                            Alpha1_rho1,
                            Alpha2_rho2,
                            Alpha1_rho1*u,
                            Alpha2_rho2*u,
#if (AMREX_SPACEDIM > 1)
                            Alpha1_rho1*v,
                            Alpha2_rho2*v,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha1_rho1*w,
                            Alpha2_rho2*w,
#endif
                            Alpha1_rho1_e1 + (0.5)*Alpha1_rho1*( AMREX_D_TERM(u*u,+v*v,+w*w) ),
                            Alpha2_rho2_e2 + (0.5)*Alpha2_rho2*( AMREX_D_TERM(u*u,+v*v,+w*w) ),
                            Alpha1,
                            Alpha2);
    
}

#endif
#endif