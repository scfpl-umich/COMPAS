// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Flux_Compressible5Eq_NPhase_H_
#define Conservative2Flux_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Flux(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                       amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                       Parm const& EOS_Parameter,
                       int ID_Dimension)
{
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real u;,
                  amrex::Real v;,
                  amrex::Real w );
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Primitive2State(Alpha_rho,
       AMREX_D_DECL(u, v, w),
                    p,
                    Alpha,
                    Primitive);

    AMREX_D_TERM( amrex::Real rho_u;,
                  amrex::Real rho_v;,
                  amrex::Real rho_w );
    amrex::Real rho_E;
    Conservative2State(Alpha_rho,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha,
                       Conservative);
    
    amrex::Real V;
    //x
    if (ID_Dimension == 0){
        V = u;
    }
#if (AMREX_SPACEDIM > 1)
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

    // if (ID_Dimension==1){
    //     p = 0.0;
    //     rho_v = 0.0;
    // }

    amrex::Array<amrex::Real,NPHASE> Flux_Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Flux_Alpha;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Flux_Alpha_rho[i_Phase] = Alpha_rho[i_Phase]*V;
        Flux_Alpha[i_Phase] = Alpha[i_Phase]*V;
    }
    Component2Vector(Flux,
                     Flux_Alpha_rho,
        AMREX_D_DECL(rho_u*V, rho_v*V, rho_w*V),
                     (rho_E + p)*V,
                     Flux_Alpha);
    Flux[NPHASE+ID_Dimension] += p;
}

#endif
#endif