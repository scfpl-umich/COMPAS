// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_SurfaceIntegralNC_Compressible6Eq_H_
#define FVM_SurfaceIntegralNC_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS == SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_SurfaceIntegral_NC(int i, int j, int k,
                            amrex::Array4<amrex::Real> const& dof_out,
                            amrex::Array4<amrex::Real const> const& dof_in,
               AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                            amrex::Array4<amrex::Real const> const& flxy,
                            amrex::Array4<amrex::Real const> const& flxz),
               AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxxNC,
                            amrex::Array4<amrex::Real const> const& flxyNC,
                            amrex::Array4<amrex::Real const> const& flxzNC),
               AMREX_D_DECL(amrex::Real dx,
                            amrex::Real dy,
                            amrex::Real dz),
                            Parm const& EOS_Parameter,
                            Physics_Parameter const& Physics_Parm)
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

    amrex::Real Div_Vel;
    FVM_Divergence_K(i, j, k, 0,
                     Div_Vel,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_Vel_Alpha1_p1;
    FVM_Divergence_K(i, j, k, 1,
                     Div_Vel_Alpha1_p1,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_Vel_Alpha2_p2;
    FVM_Divergence_K(i, j, k, 2,
                     Div_Vel_Alpha2_p2,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Real Div_Alpha1_Vel;
    FVM_Divergence_K(i, j, k, INDEX_VolumeFraction1,
                     Div_Alpha1_Vel,
        AMREX_D_DECL(flxx, flxy, flxz),
        AMREX_D_DECL(dx  , dy  , dz  ));

    amrex::Real rho = Alpha1_rho1 + Alpha2_rho2;
    amrex::Real Vel_Dot_Grad_Alpha1_p1 = Div_Vel_Alpha1_p1 - Alpha1_p1*Div_Vel;
    amrex::Real Vel_Dot_Grad_Alpha2_p2 = Div_Vel_Alpha2_p2 - Alpha2_p2*Div_Vel;

    dof_out(i,j,k,INDEX_Energy1)         +=  (Alpha2_rho2/rho*Vel_Dot_Grad_Alpha1_p1 - Alpha1_rho1/rho*Vel_Dot_Grad_Alpha2_p2);
    dof_out(i,j,k,INDEX_Energy2)         += -(Alpha2_rho2/rho*Vel_Dot_Grad_Alpha1_p1 - Alpha1_rho1/rho*Vel_Dot_Grad_Alpha2_p2);
    dof_out(i,j,k,INDEX_VolumeFraction1) += Alpha1*Div_Vel;

}


#endif
#endif