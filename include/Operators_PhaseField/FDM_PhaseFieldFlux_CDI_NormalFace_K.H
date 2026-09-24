// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_CDI_NormalFace_K_H_
#define FDM_PhaseFieldFlux_CDI_NormalFace_K_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_CDI_NormalFace_K(int i, int j, int k, int n,
                                         amrex::Real& PF_Flux,
                                         amrex::Array4<amrex::Real const> const& PF,
                                         amrex::Array4<amrex::Real const> const& Mob,
                                         amrex::Real const& Eta,
                            AMREX_D_DECL(amrex::Real const& dx,
                                         amrex::Real const& dy,
                                         amrex::Real const& dz),
                                         int ID_Dimension)
{
    amrex::Real PF_c = PF(i,j,k,n);
    amrex::Real PF_xm = PF(i-1,j,k,n);
#if (AMREX_SPACEDIM > 1)
    amrex::Real PF_xp = PF(i+1,j,k,n);
    amrex::Real PF_ym = PF(i,j-1,k,n);
    amrex::Real PF_yp = PF(i,j+1,k,n);
    amrex::Real PF_xm_ym = PF(i-1,j-1,k,n);
    amrex::Real PF_xm_yp = PF(i-1,j+1,k,n);
    amrex::Real PF_xp_ym = PF(i+1,j-1,k,n);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real PF_zm = PF(i,j,k-1,n);
    amrex::Real PF_zp = PF(i,j,k+1,n);
    amrex::Real PF_xm_zm = PF(i-1,j,k-1,n);
    amrex::Real PF_xm_zp = PF(i-1,j,k+1,n);
    amrex::Real PF_xp_zm = PF(i+1,j,k-1,n);
    amrex::Real PF_ym_zm = PF(i,j-1,k-1,n);
    amrex::Real PF_ym_zp = PF(i,j-1,k+1,n);
    amrex::Real PF_yp_zm = PF(i,j+1,k-1,n);
#endif
    
    amrex::Real Mob_Face;
    Interpolation2Face_Linear_K(i, j, k, 0,
                                Mob_Face,
                                Mob,
                                ID_Dimension);
    
    FDM_PhaseFieldFlux_CDI_NormalFace_D(PF_Flux,
                                        PF_c,
                                        PF_xm,
#if (AMREX_SPACEDIM > 1)
                                        PF_xp,
                                        PF_ym, PF_yp, PF_xm_ym, PF_xm_yp, PF_xp_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                                        PF_zm, PF_zp, PF_xm_zm, PF_xm_zp, PF_xp_zm, PF_ym_zm, PF_ym_zp, PF_yp_zm,
#endif
                                        Mob_Face, Eta,
                           AMREX_D_DECL(dx, dy, dz),
                                        ID_Dimension);
}

#endif