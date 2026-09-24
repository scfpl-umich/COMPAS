// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2ConservativePhase_Compressible5Eq_H_
#define Conservative2ConservativePhase_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2ConservativePhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2>& ConservativePhase,
                                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(amrex::Real u;, 
                 amrex::Real v;, 
                 amrex::Real w );
    amrex::Real p;
    amrex::Real Alpha1;
    amrex::Real Alpha2;
    Primitive2State(Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha1,
                    Primitive);
    Alpha2 = 1.0 - Alpha1;
    
    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Pressure to rho_e
    int ID_Phase;
    ID_Phase = 1;
    amrex::Real Alpha1_rho1_e1;
    EquationOfState_Phase(Alpha1_rho1_e1,
                          p,
                          Alpha1,
                          Alpha1_rho1,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
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
    ID_Phase = 2;
    amrex::Real Alpha2_rho2_e2;
    EquationOfState_Phase(Alpha2_rho2_e2,
                          p,
                          Alpha2,
                          Alpha2_rho2,
                          EOS_Parameter,
                          ID_Phase,
                          PressureBase,
                          OutputInternalEnergy);
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
    amrex::Real rho_e = Alpha1_rho1_e1
                      + Alpha2_rho2_e2;
    if (Alpha1_rho1_e1 < Alpha1_rho1_e1_min){
        Alpha1_rho1_e1 = Alpha1_rho1_e1_min;
        Alpha2_rho2_e2 = rho_e - Alpha1_rho1_e1;
    }
    else if (Alpha2_rho2_e2 < Alpha2_rho2_e2_min){
        Alpha2_rho2_e2 = Alpha2_rho2_e2_min;
        Alpha1_rho1_e1 = rho_e - Alpha2_rho2_e2;
    }
    
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