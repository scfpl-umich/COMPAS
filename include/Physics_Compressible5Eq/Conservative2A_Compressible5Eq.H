// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2A_Compressible5Eq_H_
#define Conservative2A_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2A(amrex::Real& A,
                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative,
                    Parm const& EOS_Parameter,
                    int ID_K)
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
    
    amrex::Array<amrex::Real,2> K;
    NonConservative(K,
                    p,
                    Alpha1,
                    EOS_Parameter,
                    ID_K);
    amrex::Real K1 = K[0];
    amrex::Real K2 = K[1];
    
    amrex::Real AP1 = EOS_Parameter.A1;
    amrex::Real AP2 = EOS_Parameter.A2;
    amrex::Real BP1 = EOS_Parameter.B1;
    amrex::Real BP2 = EOS_Parameter.B2;
    amrex::Real p_min1 = EOS_Parameter.p_min[0];
    amrex::Real p_min2 = EOS_Parameter.p_min[1];
    A = Alpha1*(1.0 + K1)*( (AP1 + 1.0)*p_min1 + BP1 )
      + Alpha2*(1.0 + K2)*( (AP2 + 1.0)*p_min2 + BP2 );
    
}

#endif
#endif