// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Pmin_Compressible5Eq_NPhase_H_
#define Conservative2Pmin_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Pmin(amrex::Real& Pmin,
                       amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                       Parm const& EOS_Parameter,
                       int ID_K)
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
    
    amrex::Array<amrex::Real,NPHASE> K;
    NonConservative(K,
                    p,
                    Alpha,
                    EOS_Parameter,
                    ID_K);
    
    amrex::Real p_mini;
    Pmin = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        p_mini = EOS_Parameter.p_min[i_Phase];
        Pmin += Alpha[i_Phase]*(1.0 + K[i_Phase])*p_mini;
    }
    
}

#endif
#endif