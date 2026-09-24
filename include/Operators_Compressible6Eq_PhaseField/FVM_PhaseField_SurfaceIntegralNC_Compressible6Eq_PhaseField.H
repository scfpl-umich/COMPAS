// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_PhaseField_SurfaceIntegralNC_Compressible6Eq_PhaseField_H_
#define FVM_PhaseField_SurfaceIntegralNC_Compressible6Eq_PhaseField_H_

#include<Macros.H>
#if (PHYSICS == SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseField_FVM_SurfaceIntegral_NC(int i, int j, int k,
                                       amrex::Array4<amrex::Real> const& dof_out,
                                       amrex::Array4<amrex::Real const> const& dof_in,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxxNC,
                                       amrex::Array4<amrex::Real const> const& flxyNC,
                                       amrex::Array4<amrex::Real const> const& flxzNC),
                          AMREX_D_DECL(amrex::Real dx,
                                       amrex::Real dy,
                                       amrex::Real dz),
                                       Parm const & EOS_Parameter)
{
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Conservative;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = dof_in(i,j,k,iState);
    }

    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);

    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM( amrex::Real u;,
                  amrex::Real v;,
                  amrex::Real w );
    amrex::Real Alpha1_p1;
    amrex::Real Alpha2_p2;
    amrex::Real Alpha1;
    Primitive2State(Alpha1_rho1,
                    Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    Alpha1_p1,
                    Alpha2_p2,
                    Alpha1,
                    Primitive);

    amrex::Real rho = Alpha1_rho1 + Alpha2_rho2;
    amrex::Real EK = 0.5*(AMREX_D_TERM(u*u, + v*v, + w*w));

    amrex::Real Div_J1_rho1;
    FVM_Divergence_K(i, j, k, 0,
                     Div_J1_rho1,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_J2_rho2;
    FVM_Divergence_K(i, j, k, 1,
                     Div_J2_rho2,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_J1_rho1_EK1;
    FVM_Divergence_K(i, j, k, 2,
                     Div_J1_rho1_EK1,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_J2_rho2_EK2;
    FVM_Divergence_K(i, j, k, 3,
                     Div_J2_rho2_EK2,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Vel_dot_J1_rho1_dot_Grad_Vel = Div_J1_rho1_EK1 - EK*Div_J1_rho1;
    amrex::Real Vel_dot_J2_rho2_dot_Grad_Vel = Div_J2_rho2_EK2 - EK*Div_J2_rho2;

    dof_out(i,j,k,INDEX_Energy1) +=  (Alpha1_rho1/rho*Vel_dot_J2_rho2_dot_Grad_Vel - Alpha2_rho2/rho*Vel_dot_J1_rho1_dot_Grad_Vel);
    dof_out(i,j,k,INDEX_Energy2) += -(Alpha1_rho1/rho*Vel_dot_J2_rho2_dot_Grad_Vel - Alpha2_rho2/rho*Vel_dot_J1_rho1_dot_Grad_Vel);

}


#endif
#endif