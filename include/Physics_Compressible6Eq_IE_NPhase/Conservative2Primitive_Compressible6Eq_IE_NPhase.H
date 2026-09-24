// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Primitive_Compressible6Eq_IE_NPhase_H_
#define Conservative2Primitive_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Primitive(amrex::Array<amrex::Real,NSTATE>& Primitive,
                            amrex::Array<amrex::Real,NSTATE> const& Conservative,
                            Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;
    AMREX_D_TERM( amrex::Real rho_u;, 
                  amrex::Real rho_v;, 
                  amrex::Real rho_w);
    amrex::Array<amrex::Real, NPHASE> Alpha_Ps;

    amrex::Real rho_E;
    amrex::Real rho = 0.0;
    amrex::Real P = 0.0;

    Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E, Conservative);


    for(int i =0; i<NPHASE; ++i){
        rho += Alpha_rho[i];
    }
    
    AMREX_D_TERM( amrex::Real u = rho_u/rho;,
                  amrex::Real v = rho_v/rho;,
                  amrex::Real w = rho_w/rho );
    
    // int ID_PT = 1;//Pressure
    // int ID_Inverse = 1;//From rho_e to Pressure

    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = false;//From rho_e to Pressure


    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        
        EquationOfState_Phase(Alpha_Ps[i_Phase],
                    Alpha_rho_e[i_Phase],
                    Alpha[i_Phase],
                    Alpha_rho[i_Phase],
                    EOS_Parameter,
                    i_Phase,
                    PressureBase,
                    OutputInternalEnergy);

        P += Alpha_Ps[i_Phase];

    }

    
    State2Primitive(Primitive, Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), P);
    
}

#endif
#endif