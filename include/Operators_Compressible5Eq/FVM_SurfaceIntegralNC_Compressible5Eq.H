// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_SurfaceIntegralNC_Compressible5Eq_H_
#define FVM_SurfaceIntegralNC_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS == FIVEEQS)

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

    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = dof_in(i,j,k,iState);
    }
    amrex::Real Alpha1 = Conservative[INDEX_VolumeFraction1];
    amrex::Array<amrex::Real,2> K;
    Conservative2K(K,
                   Conservative, EOS_Parameter,
                   ID_K);
    amrex::Real K1 = K[0];

    amrex::Real Div_Vel;
    FVM_Divergence_K(i, j, k, 0,
                     Div_Vel,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));

    dof_out(i,j,k,INDEX_VolumeFraction1) += Alpha1*(1.0 + K1)*Div_Vel;
}


#endif
#endif