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

#if (DIFFUSION == true)
    // Viscous energy of the phases. The conservative flux puts all of [div(u.tau)]_h into
    // phase 1; phase 2 takes Y_2 [div(u.tau)]_h + alpha_2 Phi_2 - Y_2 Phi from it, with
    // Phi_k = tau_k : grad u and Phi = alpha_1 Phi_1 + alpha_2 Phi_2. Each phase then gets
    // Y_k [div(u.tau)]_h + alpha_k Phi_k - Y_k Phi: the sum is conservative, an absent phase
    // (alpha_k = Y_k = 0) gets 0 and a uniform velocity gives 0.
    amrex::Real Div_Work;
    FVM_Divergence_K(i, j, k, INDEX_NC_ViscousWork,
                     Div_Work,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    amrex::Array<amrex::Real,2> Phi_Phase;
    FDM_Conservative2ViscousDissipation_K(i, j, k,
                                          Phi_Phase,
                                          dof_in,
                             AMREX_D_DECL(dx, dy, dz),
                                          EOS_Parameter);

    amrex::Real Alpha2   = 1.0 - Alpha1;
    amrex::Real Y2       = Alpha2_rho2/rho;
    amrex::Real Phi      = Alpha1*Phi_Phase[0] + Alpha2*Phi_Phase[1];
    amrex::Real Transfer = Y2*Div_Work + (Alpha2*Phi_Phase[1] - Y2*Phi);

    dof_out(i,j,k,INDEX_Energy1) -= Transfer;
    dof_out(i,j,k,INDEX_Energy2) += Transfer;
#endif

}


#endif
#endif