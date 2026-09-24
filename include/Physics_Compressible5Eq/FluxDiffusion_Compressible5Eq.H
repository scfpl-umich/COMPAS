// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxDiffusion_Compressible5Eq_H_
#define FluxDiffusion_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxDiffusion(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                   amrex::Real& dmax,
                   amrex::Array<amrex::Real,AMREX_SPACEDIM> const& Vel,
                   amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Grad_Vel,
                   amrex::Real const& Grad_T,
                   amrex::Real const& mu,
                   amrex::Real const& mub,
                   amrex::Real const& Kappa,
                   int ID_Dimension)
{
    amrex::Real tau_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Real tau_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real tau_w;
#endif
    amrex::Real tau_b;
    amrex::Real Work;
    amrex::Real qC;
    
    tau_u = mu*(Grad_Vel[0][ID_Dimension] + Grad_Vel[ID_Dimension][0]);
    Work = Vel[0]*tau_u;
#if (AMREX_SPACEDIM > 1)
    tau_v = mu*(Grad_Vel[1][ID_Dimension] + Grad_Vel[ID_Dimension][1]);
    Work += Vel[1]*tau_v;
#endif
#if (AMREX_SPACEDIM > 2)
    tau_w = mu*(Grad_Vel[2][ID_Dimension] + Grad_Vel[ID_Dimension][2]);
    Work += Vel[2]*tau_w;
#endif
    amrex::Real Div_Vel = AMREX_D_TERM( Grad_Vel[0][0],
                                      + Grad_Vel[1][1],
                                      + Grad_Vel[2][2] );
    tau_b = (mub - 2.0/3.0*mu)*Div_Vel;
    Work += Vel[ID_Dimension]*tau_b;
    qC = Kappa*Grad_T;
    
    //Diffusion flux
    Component2Vector(Flux,
                     0.0,
                     0.0,
        AMREX_D_DECL(tau_u, tau_v, tau_w),
                     Work + qC,
                     0.0);
    Flux[2+ID_Dimension] += tau_b;
    
    dmax = mu + (mub - 2.0/3.0*mu) + Kappa;
    
}

#endif
#endif