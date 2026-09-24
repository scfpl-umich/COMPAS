// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_K_Compressible5Eq_H_
#define FDM_Conservative2FluxDiffusion_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_K(int i, int j, int k,
                                      amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                                      amrex::Real& dmax,
                                      amrex::Array4<amrex::Real const> const& Conservative,
                                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                         AMREX_D_DECL(amrex::Real const& dx,
                                      amrex::Real const& dy,
                                      amrex::Real const& dz),
                                      amrex::Real const& time,
                                      Parm const& EOS_Parameter,
                                      int ID_Dimension)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_c;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xm;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_ym;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_yp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xm_ym;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xm_yp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xp_ym;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_zm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_zp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xm_zm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xm_zp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_xp_zm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_ym_zm;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_ym_zp;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_yp_zm;
#endif
    for (int i_State=0; i_State<4+AMREX_SPACEDIM; i_State++){
        Conservative_c[i_State] = Conservative(i,j,k,i_State);
        Conservative_xm[i_State] = Conservative(i-1,j,k,i_State);
#if (AMREX_SPACEDIM > 1)
        Conservative_xp[i_State] = Conservative(i+1,j,k,i_State);
        Conservative_ym[i_State] = Conservative(i,j-1,k,i_State);
        Conservative_yp[i_State] = Conservative(i,j+1,k,i_State);
        Conservative_xm_ym[i_State] = Conservative(i-1,j-1,k,i_State);
        Conservative_xm_yp[i_State] = Conservative(i-1,j+1,k,i_State);
        Conservative_xp_ym[i_State] = Conservative(i+1,j-1,k,i_State);
#endif
#if (AMREX_SPACEDIM > 2)
        Conservative_zm[i_State] = Conservative(i,j,k-1,i_State);
        Conservative_zp[i_State] = Conservative(i,j,k+1,i_State);
        Conservative_xm_zm[i_State] = Conservative(i-1,j,k-1,i_State);
        Conservative_xm_zp[i_State] = Conservative(i-1,j,k+1,i_State);
        Conservative_xp_zm[i_State] = Conservative(i+1,j,k-1,i_State);
        Conservative_ym_zm[i_State] = Conservative(i,j-1,k-1,i_State);
        Conservative_ym_zp[i_State] = Conservative(i,j-1,k+1,i_State);
        Conservative_yp_zm[i_State] = Conservative(i,j+1,k-1,i_State);
#endif
    }
    amrex::Real mu;
    amrex::Real mub;
    Conservative2Viscosity_Face_K(i, j, k,
                                  mu,
                                  mub,
                                  Conservative,
                                  prob_lo,
                     AMREX_D_DECL(dx, dy, dz),
                                  time,
                                  EOS_Parameter,
                                  ID_Dimension);
    amrex::Real Kappa;
    Conservative2Conductivity_Face_K(i, j, k,
                                     Kappa,
                                     Conservative,
                                     prob_lo,
                        AMREX_D_DECL(dx, dy, dz),
                                     time,
                                     EOS_Parameter,
                                     ID_Dimension);

    
    FDM_Conservative2FluxDiffusion_D(Flux,
                                     dmax,
                                     Conservative_c,
                                     Conservative_xm,
#if (AMREX_SPACEDIM > 1)
                                     Conservative_xp, Conservative_ym, Conservative_yp,
                                     Conservative_xm_ym, Conservative_xm_yp, Conservative_xp_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                                     Conservative_zm, Conservative_zp,
                                     Conservative_xm_zm, Conservative_xm_zp, Conservative_xp_zm,
                                     Conservative_ym_zm, Conservative_ym_zp, Conservative_yp_zm,
#endif
                                     mu, mub, Kappa,
                        AMREX_D_DECL(dx, dy, dz),
                                     EOS_Parameter,
                                     ID_Dimension);
    
}

#endif
#endif

