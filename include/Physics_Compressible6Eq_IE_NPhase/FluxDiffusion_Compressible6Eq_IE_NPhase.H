// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxDiffusion_Compressible6Eq_IE_NPhase_H_
#define FluxDiffusion_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxDiffusion(amrex::Array<amrex::Real,NSTATE>& Flux,
                   amrex::Real& dmax,
                   amrex::Array<amrex::Real,AMREX_SPACEDIM>& tau,
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
    tau[0] = tau_u;
    // amrex::Print()<<"mu: "<<mu<<" mub: "<<mub<<" Kappa: "<<Kappa<<"\n";
    // amrex::Print()<<"Vel u: "<<Vel[0]<<" v: "<<Vel[1]<<"\n";
    // amrex::Print()<<"GradV [0][0]: "<<Grad_Vel[0][0]<<" GradV [0][1]: "<<Grad_Vel[0][1]<<"\n";
    // amrex::Print()<<"GradV [1][0]: "<<Grad_Vel[1][0]<<" GradV [1][1]: "<<Grad_Vel[1][1]<<"\n";
#if (AMREX_SPACEDIM > 1)
    tau_v = mu*(Grad_Vel[1][ID_Dimension] + Grad_Vel[ID_Dimension][1]);
    Work += Vel[1]*tau_v;
    tau[1] = tau_v;
#endif
#if (AMREX_SPACEDIM > 2)
    tau_w = mu*(Grad_Vel[2][ID_Dimension] + Grad_Vel[ID_Dimension][2]);
    Work += Vel[2]*tau_w;
    tau[2] = tau_w;
#endif
    amrex::Real Div_Vel = AMREX_D_TERM( Grad_Vel[0][0],
                                      + Grad_Vel[1][1],
                                      + Grad_Vel[2][2] );
    tau_b = (mub - 2.0/3.0*mu)*Div_Vel;
    Work += Vel[ID_Dimension]*tau_b;
    qC = Kappa*Grad_T;
    
    //Diffusion flux
    for(int i_state=0; i_state < NSTATE; ++i_state){
        Flux[i_state] = 0.;
    }

    Flux[ru] = tau_u;
    
    // Flux[2] = tau_u;
#if (AMREX_SPACEDIM > 1)
    // Flux[3] = tau_v;
    Flux[rv] = tau_v;
#endif
#if (AMREX_SPACEDIM > 2)
    // Flux[4] = tau_w;
    Flux[rw] = tau_w;
#endif
    Flux[rE] = Work + qC;

    Flux[ru+ID_Dimension] += tau_b;

    // amrex::Print()<<"tau_u "<<tau_u<<" tau_v: "<<tau_v<<" Work "<<Work<<" qC "<<qC<<"\n";
    
    dmax = mu + (mub - 2.0/3.0*mu) + Kappa;
    
}

#endif
#endif