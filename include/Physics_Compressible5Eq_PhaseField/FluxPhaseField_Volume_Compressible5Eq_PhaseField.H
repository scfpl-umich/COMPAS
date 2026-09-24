// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Volume_Compressible5Eq_PhaseField_H_
#define FluxPhaseField_Volume_Compressible5Eq_PhaseField_H_


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
PhaseFieldVolumeFlux (
amrex::Array<amrex::Array<amrex::Real,4+AMREX_SPACEDIM>, AMREX_SPACEDIM> & Flux,
amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const & Conservative,
amrex::Array<amrex::Real, AMREX_SPACEDIM> const & GradQ,
Parm const& EOS_Parameter)
{

    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(
    amrex::Real rho_u;,
    amrex::Real rho_v;,
    amrex::Real rho_w;);
    amrex::Real rho_E;
    amrex::Real Alpha1;

    AMREX_D_TERM(
    amrex::Real u;,
    amrex::Real v;,
    amrex::Real w;);
    amrex::Real p;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Primitive;

    Conservative2Primitive(Primitive, Conservative, EOS_Parameter);
    Conservative2State(Alpha1_rho1,Alpha2_rho2,AMREX_D_DECL(rho_u,rho_v,rho_w),rho_E,Alpha1,Conservative);
    Primitive2State(Alpha1_rho1,Alpha2_rho2,AMREX_D_DECL(u,v,w),p,Alpha1,Primitive);

    amrex::Array<amrex::Real, AMREX_SPACEDIM> Grad_Q1;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> Grad_Q2;
    for (int d = 0; d < AMREX_SPACEDIM; d++){
        Grad_Q1[d] = GradQ[d]; Grad_Q2[d] = -GradQ[d];
    }

    amrex::Real Alpha2 = 1.0 - Alpha1;

    for (int d = 0; d < AMREX_SPACEDIM; d++){
        

        amrex::Real J1 = Alpha1*Alpha2*Grad_Q1[d];
        amrex::Real J2 = Alpha2*Alpha1*Grad_Q2[d];
        amrex::Real J1_rho1 = Alpha1_rho1*Alpha2*Grad_Q1[d];
        amrex::Real J2_rho2 = Alpha2_rho2*Alpha1*Grad_Q2[d];


        //Momentum
        amrex::Real M = J1_rho1 + J2_rho2;
        
        amrex::Real M_u = M*u;
#if (AMREX_SPACEDIM > 1)
        amrex::Real M_v = M*v;
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real M_w = M*w;
#endif

        //Energy
        amrex::Real eK = 0.5*( AMREX_D_TERM(u*u,+v*v,+w*w) );
        amrex::Real M_eK = M*eK;
        amrex::Real J1_rho1_e1 = EOS_Parameter.A1*J1*p + EOS_Parameter.B1*J1 + EOS_Parameter.D1*J1_rho1;
        amrex::Real J2_rho2_e2 = EOS_Parameter.A2*J2*p + EOS_Parameter.B2*J2 + EOS_Parameter.D2*J2_rho2;

        
        //Flux_PhaseField
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux_d = Flux[d];
        Component2Vector(Flux_d,
                         J1_rho1, J2_rho2,
            AMREX_D_DECL(M_u, M_v, M_w),
                         M_eK + J1_rho1_e1 + J2_rho2_e2,
                         J1);
        
        //Flux[d][0] = J1_rho1;
        //Flux[d][1] = J2_rho2;
        //Flux[d][2] = M_u;
//#if (AMREX_SPACEDIM > 1)
        //Flux[d][3] = M_v;
//#endif
//#if (AMREX_SPACEDIM > 2)
        //Flux[d][4] = M_w;
//#endif
        //Flux[d][2+AMREX_SPACEDIM] = M_eK + J1_rho1_e1 + J2_rho2_e2;
        //Flux[d][3+AMREX_SPACEDIM] = J1;
    }
    
}

#endif