// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2A_Compressible5Eq_NPhase_H_
#define Conservative2A_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2A(amrex::Real& A,
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
                    Alpha_rho,
                    EOS_Parameter,
                    ID_K);
    
    amrex::Real EOS_APi;
    amrex::Real EOS_BPi;
    amrex::Real p_mini;
    A = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        EOS_APi = EOS_Parameter.A[i_Phase];
        EOS_BPi = EOS_Parameter.B[i_Phase];
        p_mini  = EOS_Parameter.p_min[i_Phase];
        // Kapila coefficient alpha_k (1 + K_k): phi_k (1 + K_k) for NASG when K is on (see NonConservative)
        amrex::Real Phi = (ID_K == 1) ? NASG_Phi(Alpha[i_Phase], Alpha_rho[i_Phase], EOS_Parameter.C[i_Phase]) : Alpha[i_Phase];
        A += Phi*(1.0 + K[i_Phase])*( (EOS_APi + 1.0)*p_mini + EOS_BPi );
    }
    
}

#endif
#endif