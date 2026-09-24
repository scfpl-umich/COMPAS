// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Lambda_HLL_Compressible5Eq_NPhase_H_
#define Conservative2Lambda_HLL_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Lambda_HLL(amrex::Real& Lambda_p,
                             amrex::Real& Lambda_m,
                             amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                             amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative0,
                             Parm const& EOS_Parameter,
                             int ID_Dimension,
                             int ID_K,
                             int ID_Lambda)
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
               Alpha,
               Alpha_rho,
               EOS_Parameter,
               ID_K);
    
    amrex::Real C = c;
    if ( ID_Lambda == 0 || std::isnan(c) ){
        amrex::Real rho = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            rho += Alpha_rho[i_Phase];
        }
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