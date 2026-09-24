// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef AnalyticDiffusiveFlux_Compressible5Eq_H_
#define AnalyticDiffusiveFlux_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

#include<Conservative2Velocity_Compressible5Eq.H>
#include<ConservativeGradients2PrimitiveGradients_Compressible5Eq.H>
#include<FluxDiffusion_Compressible5Eq.H>

template<typename TF, typename TD, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AnalyticDiffusiveFlux (TF & Ghat,
					   AMREX_D_DECL(
					   TD const & dUdXi,
					   TD const & dUdEta,
					   TD const & dUdZeta),
					   TU const & Conservative,
					   amrex::Real & dmax,
					   int ID_Dim,
					   Parm const & EOS_Parameter)
{


    AMREX_D_TERM( amrex::Real u;,
                  amrex::Real v;,
                  amrex::Real w);
    Conservative2Velocity(AMREX_D_DECL(u, v, w),Conservative);

    amrex::Array<amrex::Array<amrex::Real, NSTATE>, AMREX_SPACEDIM> U_Grad = {AMREX_D_DECL(dUdXi, dUdEta, dUdZeta)};

    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM>,AMREX_SPACEDIM> Vel_Grad;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> Vel = {AMREX_D_DECL(u,v,w)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> T_Grad;

    ConservativeGradients2PrimitiveGradients(Vel_Grad,T_Grad,U_Grad,Conservative,EOS_Parameter);
                                              
    amrex::Real mu;
    amrex::Real mub;
    Conservative2Viscosity(mu,mub,Conservative,EOS_Parameter);
    amrex::Real Kappa;
    Conservative2Conductivity(Kappa,Conservative,EOS_Parameter);

	FluxDiffusion(Ghat,dmax,Vel,Vel_Grad,T_Grad[ID_Dim],mu,mub,Kappa,ID_Dim);


}


template<typename TF, typename TD, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AnalyticDiffusiveFlux (AMREX_D_DECL(
					    TF & FluxX,
					    TF & FluxY,
					    TF & FluxZ),
						AMREX_D_DECL(
					    TD const & dUdXi,
					    TD const & dUdEta,
					    TD const & dUdZeta),
						TU const & Conservative,
						Parm const & parm)
{
	amrex::Real c;

	AnalyticDiffusiveFlux(FluxX,AMREX_D_DECL(dUdXi,dUdEta,dUdZeta),Conservative, c, 0, parm);
#if (AMREX_SPACEDIM > 1)
	AnalyticDiffusiveFlux(FluxY,AMREX_D_DECL(dUdXi,dUdEta,dUdZeta),Conservative, c, 1, parm);
#endif
#if (AMREX_SPACEDIM > 2)
	AnalyticDiffusiveFlux(FluxZ,AMREX_D_DECL(dUdXi,dUdEta,dUdZeta),Conservative, c, 2, parm);
#endif


}

#endif
#endif