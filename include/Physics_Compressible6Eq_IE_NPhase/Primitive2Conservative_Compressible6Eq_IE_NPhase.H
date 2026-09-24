// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Primitive2Conservative_Compressible6Eq_IE_NPhase_H_
#define Primitive2Conservative_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Primitive2Conservative(amrex::Array<amrex::Real, NSTATE>& Conservative,
                            amrex::Array<amrex::Real, NSTATE> const& Primitive,
                            Parm const& EOS_Parameter)
{
   amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;

    AMREX_D_TERM( amrex::Real u;, 
                  amrex::Real v;, 
                  amrex::Real w);

    amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
    amrex::Real P;

    Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), P, Primitive);
   
    amrex::Real rho = 0.0;

    // int ID_PT = 1;//Pressure
    // int ID_Inverse=0;//From Pressure to rho_e

    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From rho_e to Pressure

    amrex::Real rho_e = 0.0;

    for(int i =0; i<NPHASE; ++i){
        rho += Alpha_rho[i];

        EquationOfState_Phase(Alpha_rho_e[i],
                    Alpha_Ps[i],
                    Alpha[i],
                    Alpha_rho[i],
                    EOS_Parameter,
                    i,
                    PressureBase,
                    OutputInternalEnergy);

        rho_e += Alpha_rho_e[i];
    }

    
    AMREX_D_TERM( amrex::Real rho_u = rho*u;,
                  amrex::Real rho_v = rho*v;,
                  amrex::Real rho_w = rho*w );
    

   
    amrex::Real rho_E;
   
    rho_E = rho_e + (0.5)*rho*( AMREX_D_TERM( u*u, +v*v, +w*w ) );
   
    
    State2Conservative(Conservative, Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E);
    
}

#endif
#endif