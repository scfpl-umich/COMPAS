// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Velocity_Compressible5Eq_H_
#define Conservative2Velocity_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Velocity(AMREX_D_DECL(amrex::Real& u,
                                        amrex::Real& v,
                                        amrex::Real& w),
                           amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative)
{
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM( amrex::Real rho_u;, 
                  amrex::Real rho_v;, 
                  amrex::Real rho_w );
    amrex::Real rho_E;
    amrex::Real Alpha1;
    Conservative2State(Alpha1_rho1,
                       Alpha2_rho2,
          AMREX_D_DECL(rho_u, rho_v, rho_w),
                       rho_E,
                       Alpha1,
                       Conservative);
    
    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;
    AMREX_D_TERM( u = rho_u/rho;,
                  v = rho_v/rho;,
                  w = rho_w/rho );
    
}

#endif
#endif