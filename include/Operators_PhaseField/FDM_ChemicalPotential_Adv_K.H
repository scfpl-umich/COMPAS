// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_ChemicalPotential_Adv_K_H_
#define FDM_ChemicalPotential_Adv_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_ChemicalPotential_Adv_K(int i, int j, int k, int n,
                                 amrex::Real& ChemicalPotential,
                                 amrex::Array4<amrex::Real const> const& PF,
                                 amrex::Array4<amrex::Real const> const& Normal,
                                 amrex::Real const& Lambda,
                                 amrex::Real const& Eta,
                    AMREX_D_DECL(amrex::Real const& dx,
                                 amrex::Real const& dy,
                                 amrex::Real const& dz),
                                 int ID_Potential = 0)
{
    amrex::Real PF_c = PF(i,j,k,n);
    amrex::Real PF_xm = PF(i-1,j,k,n);
    amrex::Real PF_xp = PF(i+1,j,k,n);
    amrex::Real Normalx_c = Normal(i,j,k,0);
    amrex::Real Normalx_xm = Normal(i-1,j,k,0);
    amrex::Real Normalx_xp = Normal(i+1,j,k,0);
#if (AMREX_SPACEDIM > 1)
    amrex::Real PF_ym = PF(i,j-1,k,n);
    amrex::Real PF_yp = PF(i,j+1,k,n);
    amrex::Real Normaly_c = Normal(i,j,k,1);
    amrex::Real Normaly_ym = Normal(i,j-1,k,1);
    amrex::Real Normaly_yp = Normal(i,j+1,k,1);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real PF_zm = PF(i,j,k-1,n);
    amrex::Real PF_zp = PF(i,j,k+1,n);
    amrex::Real Normalz_c = Normal(i,j,k,2);
    amrex::Real Normalz_zm = Normal(i,j,k-1,2);
    amrex::Real Normalz_zp = Normal(i,j,k+1,2);
#endif
    FDM_ChemicalPotential_Adv_D(ChemicalPotential,
                                PF_c,
                                PF_xm, PF_xp,
                                Normalx_c, Normalx_xm, Normalx_xp,
#if (AMREX_SPACEDIM > 1)
                                PF_ym, PF_yp,
                                Normaly_c, Normaly_ym, Normaly_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                PF_zm, PF_zp,
                                Normalz_c, Normalz_zm, Normalz_zp,
#endif
                                Lambda, Eta,
                   AMREX_D_DECL(dx, dy, dz),
                                ID_Potential);
    
}

#endif