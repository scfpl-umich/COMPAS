// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2CFLSpeed_HLL_Compressible5Eq_NPhase_H_
#define Conservative2CFLSpeed_HLL_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2CFLSpeed_HLL(amrex::Real& CFLSpeed,
                               amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative,
                               amrex::Real const& cmax,
                               amrex::Real const& Quadrature_Weight,
                               Parm const& EOS_Parameter,
                               int ID_K)
{
    amrex::Array<amrex::Real,NPHASE> K;
    Conservative2K(K,
                   Conservative,
                   EOS_Parameter,
                   ID_K);
    
    amrex::Real max_abs_K = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        max_abs_K = std::max( max_abs_K, std::abs(1.0 + K[i_Phase]) );
    }

    AMREX_D_TERM(amrex::Real u;, amrex::Real v;, amrex::Real w);
    Conservative2Velocity(AMREX_D_DECL(u, v, w),
                          Conservative);
    //amrex::Real V = std::sqrt(AMREX_D_TERM(u*u, + v*v, + w*w));
    amrex::Real V = cmax;
    CFLSpeed = 2.0*(cmax + Quadrature_Weight*max_abs_K*V);
    
}

#endif
#endif