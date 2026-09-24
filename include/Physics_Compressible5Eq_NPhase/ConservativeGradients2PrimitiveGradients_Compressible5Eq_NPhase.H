// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativeGradients2PrimitiveGradients_Compressible5Eq_NPhase_H_
#define ConservativeGradients2PrimitiveGradients_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


template<typename TD>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativeGradients2PrimitiveGradients(amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM>,AMREX_SPACEDIM>& Vel_Grad,
                                              amrex::Array<amrex::Real, AMREX_SPACEDIM>& T_Grad,
                                              TD U_Grad,
                                              amrex::Array<amrex::Real, NSTATE> const& Conservative,
                                              Parm const& EOS_Parameter)
{

  
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;
    Conservative2Primitive(Primitive,Conservative,EOS_Parameter);

    AMREX_D_TERM(amrex::Real u;,
                 amrex::Real v;,
                 amrex::Real w);
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    amrex::Real p;

    Primitive2State(Alpha_rho,AMREX_D_DECL(u,v,w),p,Alpha,Primitive);

    amrex::Real T;
    Conservative2Temperature(T,Conservative,EOS_Parameter);


    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,NPHASE> Alpha_Grad;
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,NPHASE> Alpha_rho_Grad;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> rho_Grad = {AMREX_D_DECL(0.0,0.0,0.0)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> rho_cv_Grad = {AMREX_D_DECL(0.0,0.0,0.0)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> rhoE_Grad;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> P_Grad;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> Num_Grad;

    amrex::Real rho = 0.0;
    amrex::Real rho_cv = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
        rho_cv += Alpha_rho[i_Phase]*EOS_Parameter.cv[i_Phase];
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Alpha_rho_Grad[i_Phase][idim] = U_Grad[idim][i_Phase];
            Alpha_Grad[i_Phase][idim]     = U_Grad[idim][i_Phase+INDEX_VolumeFraction1];
            rho_Grad[idim] += Alpha_rho_Grad[i_Phase][idim];
            rho_cv_Grad[idim] += Alpha_rho_Grad[i_Phase][idim]*EOS_Parameter.cv[i_Phase];
        }

        rhoE_Grad[idim] = U_Grad[idim][INDEX_Energy];

        AMREX_D_TERM(
        Vel_Grad[0][idim] = (U_Grad[idim][INDEX_MomentumX] - u*rho_Grad[idim])/rho;,
        Vel_Grad[1][idim] = (U_Grad[idim][INDEX_MomentumY] - v*rho_Grad[idim])/rho;,
        Vel_Grad[2][idim] = (U_Grad[idim][INDEX_MomentumZ] - w*rho_Grad[idim])/rho;
        );


        amrex::Real Num = Conservative[INDEX_Energy] - 0.5*rho*(AMREX_D_TERM(u*u,+v*v,+w*w));
        Num_Grad[idim] = rhoE_Grad[idim] 
                         - 0.5*rho_Grad[idim]*(AMREX_D_TERM(u*u,+v*v,+w*w))
                         - rho*(AMREX_D_TERM(u*Vel_Grad[0][idim],+v*Vel_Grad[1][idim],+w*Vel_Grad[2][idim]));

        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){

            Num += - Alpha[i_Phase]*EOS_Parameter.B[i_Phase]/(1. + EOS_Parameter.A[i_Phase]) 
                   + Alpha_rho[i_Phase]*EOS_Parameter.D[i_Phase];

            Num_Grad[idim] += - Alpha_Grad[i_Phase][idim]*EOS_Parameter.B[i_Phase]/(1. + EOS_Parameter.A[i_Phase]) 
                              + Alpha_rho_Grad[i_Phase][idim]*EOS_Parameter.D[i_Phase];
        }

        T_Grad[idim] = (Num_Grad[idim] - Num*rho_cv_Grad[idim]/rho_cv)/rho_cv;


    }

    
}

#endif
#endif