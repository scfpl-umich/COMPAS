// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_SurfaceIntegralNC_Compressible5Eq_NPhase_H_
#define FVM_SurfaceIntegralNC_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


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
    int ID_K = Physics_Parm.source_term;

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative;
    for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
        Conservative[i_State] = dof_in(i,j,k,i_State);
    }

    amrex::Array<amrex::Real,NPHASE> K;
    Conservative2K(K,
                   Conservative, EOS_Parameter,
                   ID_K);

    amrex::Real Div_Vel;
    FVM_Divergence_K(i, j, k, 0,
                     Div_Vel,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::Real Alphai = Conservative[INDEX_VolumeFraction1+i_Phase];
        // Kapila coefficient alpha_k (1 + K_k) = phi_k (1 + K_k) for NASG (see NonConservative)
        amrex::Real Phii = (ID_K == 1) ? NASG_Phi(Alphai, Conservative[INDEX_Mass1+i_Phase], EOS_Parameter.C[i_Phase]) : Alphai;
        dof_out(i,j,k,INDEX_VolumeFraction1+i_Phase) += Phii*(1.0 + K[i_Phase])*Div_Vel;
    }
    
}


#endif
#endif