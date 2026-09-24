// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Flux_Compressible6Eq_IE_NPhase_H_
#define Conservative2Flux_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Flux(amrex::Array<amrex::Real,NSTATE>& Flux,
                       amrex::Array<amrex::Real,NSTATE> const& Conservative,
                       Parm const& EOS_Parameter,
                       int ID_Dimension)
{
    amrex::Array<amrex::Real, NSTATE> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real u;, 
                  amrex::Real v;, 
                  amrex::Real w);
    amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
    amrex::Real p;

    Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), p, Primitive);

    AMREX_D_TERM( amrex::Real rho_u;,
                  amrex::Real rho_v;,
                  amrex::Real rho_w );

    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;
    amrex::Real rho_E;

    Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E, Conservative);

    amrex::Real V;
    //x
    if (ID_Dimension == 0){
        V = u;
    }
    #if(AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        V = v;
    }
    #endif
    #if (AMREX_SPACEDIM > 2)
    //z
    else if (ID_Dimension == 2){
        V = w;
    }
    #endif

    for (int i_state = 0; i_state < NSTATE; ++i_state){
        Flux[i_state] = Conservative[i_state]*V;
    }

    //x
    if (ID_Dimension == 0){

        Flux[ru] += p; // rho_u
    } 

    #if(AMREX_SPACEDIM > 1)
    // y
    else if (ID_Dimension == 1){

        Flux[rv] += p; // rho_v
    } 
    #endif
    #if(AMREX_SPACEDIM > 2)
    // z
    else if (ID_Dimension == 2){

        Flux[rw] += p; // rho_w
    }
    #endif

    Flux[rE] += p*V; // rho_E
}

#endif
#endif