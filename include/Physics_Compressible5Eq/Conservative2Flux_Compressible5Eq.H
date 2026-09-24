// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Flux_Compressible5Eq_H_
#define Conservative2Flux_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Flux(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                       Parm const& EOS_Parameter,
                       int ID_Dimension)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM( amrex::Real u;,
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

    AMREX_D_TERM( amrex::Real rho_u;,
                  amrex::Real rho_v;,
                  amrex::Real rho_w );
    amrex::Real rho_E;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha1,
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

    Component2Vector(Flux,
                     Alpha1_rho1*V,
                     Alpha2_rho2*V,
        AMREX_D_DECL(rho_u*V, rho_v*V, rho_w*V),
                     (rho_E + p)*V,
                     Alpha1*V);
    Flux[2+ID_Dimension] += p;
}

#endif
#endif