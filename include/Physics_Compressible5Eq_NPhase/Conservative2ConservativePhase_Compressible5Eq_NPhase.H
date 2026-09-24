// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2ConservativePhase_Compressible5Eq_NPhase_H_
#define Conservative2ConservativePhase_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2ConservativePhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& ConservativePhase,
                                    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
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
    
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e_min;
    int isNegative = 0;
    amrex::Real rho_e = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho_u[i_Phase] = Alpha_rho[i_Phase]*u;
#if (AMREX_SPACEDIM > 1)
        Alpha_rho_v[i_Phase] = Alpha_rho[i_Phase]*v;
#endif
#if (AMREX_SPACEDIM > 2)
        Alpha_rho_w[i_Phase] = Alpha_rho[i_Phase]*w;
#endif
        bool PressureBase = true;//Pressure
        bool OutputInternalEnergy = true;//From Pressure to Alpha_rho_e
        amrex::Real Alphai_rhoi_ei;
        EquationOfState_Phase(Alphai_rhoi_ei,
                              p,
                              Alpha[i_Phase],
                              Alpha_rho[i_Phase],
                              EOS_Parameter,
                              i_Phase,
                              PressureBase,
                              OutputInternalEnergy);
        Alpha_rho_e[i_Phase] = Alphai_rhoi_ei;
        rho_e += Alphai_rhoi_ei;
        
        amrex::Real p_mini = EOS_Parameter.p_min[i_Phase];
        amrex::Real Alphai_rhoi_ei_min;
        EquationOfState_Phase(Alphai_rhoi_ei_min,
                              p_mini,
                              Alpha[i_Phase],
                              Alpha_rho[i_Phase],
                              EOS_Parameter,
                              i_Phase,
                              PressureBase,
                              OutputInternalEnergy);
        Alpha_rho_e_min[i_Phase] = Alphai_rhoi_ei_min;
        
        if (Alphai_rhoi_ei < Alphai_rhoi_ei_min){
            isNegative = 1;
        }
    }
    
    amrex::Array<int,NPHASE> Phase_Indicator = {0};
    while (isNegative){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Alpha_rho_e[i_Phase] < Alpha_rho_e_min[i_Phase]){
                Alpha_rho_e[i_Phase] = Alpha_rho_e_min[i_Phase];
                Phase_Indicator[i_Phase] = 1;
            }
        }
        
        amrex::Real Sum_Alpha_rho_e = 0.0;
        amrex::Real Sum_Alpha = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Sum_Alpha_rho_e += Alpha_rho_e[i_Phase];
            if (Phase_Indicator[i_Phase] == 0){
                Sum_Alpha += Alpha[i_Phase];
            }
        }
        
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Phase_Indicator[i_Phase] == 0){
                Alpha_rho_e[i_Phase] += Alpha[i_Phase]/Sum_Alpha*(rho_e - Sum_Alpha_rho_e);
            }
        }
        
        isNegative = 0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Alpha_rho_e[i_Phase] < Alpha_rho_e_min[i_Phase]){
                isNegative = 1;
            }
        }
    }
    
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_E;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho_E[i_Phase] = (0.5)*Alpha_rho[i_Phase]*( AMREX_D_TERM(u*u,+v*v,+w*w) ) + Alpha_rho_e[i_Phase];
    }
    State2ConservativePhase(ConservativePhase,
                            Alpha_rho,
                            Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                            Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha_rho_w,
#endif
                            Alpha_rho_E,
                            Alpha);
    
}

#endif
#endif