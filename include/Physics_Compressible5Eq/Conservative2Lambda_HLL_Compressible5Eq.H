// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Lambda_HLL_Compressible5Eq_H_
#define Conservative2Lambda_HLL_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Lambda_HLL(amrex::Real& Lambda_p,
                             amrex::Real& Lambda_m,
                             amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                             amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative0,
                             Parm const& EOS_Parameter,
                             int ID_Dimension,
                             int ID_K,
                             int ID_Lambda)
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
    
    amrex::Real c;
    SoundSpeed(c,
               p,
               Alpha1,
               Alpha1_rho1,
               Alpha2_rho2,
               EOS_Parameter,
               ID_K);
    
    amrex::Real C = c;
    if ( ID_Lambda == 0 || std::isnan(c) ){
        amrex::Real rho = Alpha1_rho1
                        + Alpha2_rho2;
        amrex::Real A0;
        Conservative2A(A0,
                       Conservative0,
                       EOS_Parameter,
                       ID_K);
        amrex::Real Pmin0;
        Conservative2Pmin(Pmin0,
                          Conservative0,
                          EOS_Parameter,
                          ID_K);
        amrex::Real EnergyBound;
        Conservative2EnergyBound(EnergyBound,
                                 Conservative,
                                 EOS_Parameter);
        C = std::max( std::abs(p - Pmin0 + A0)/std::sqrt(2.0*rho*EnergyBound), c );
    }
    
    Lambda_p = V + C;
    Lambda_m = V - C;
    
}

#endif
#endif