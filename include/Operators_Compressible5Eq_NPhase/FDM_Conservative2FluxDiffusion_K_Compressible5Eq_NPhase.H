// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_K_Compressible5Eq_NPhase_H_
#define FDM_Conservative2FluxDiffusion_K_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_K(int i, int j, int k,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
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
    // Face stencil: 1 ghost cell with the edge neighbors. FaceStencil::OffsetOnFace maps the points
    // one cell above the face cell along the face normal (xp, xp_ym, xp_zm on an x face), which this
    // face does not use, to the face cell, so that no read leaves the single ghost layer of
    // FiniteVolume.Scheme = GODUNOV (see FDM_FaceVelocityGradient_D.H). The values used are unchanged.
    using namespace FaceStencil;
    amrex::Array<amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>,FaceStencil::N> Conservative_s;
    for (int s = 0; s < FaceStencil::N; ++s){
        int di, dj, dk;
        FaceStencil::OffsetOnFace(s, ID_Dimension, di, dj, dk);
        for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
            Conservative_s[s][i_State] = Conservative(i+di,j+dj,k+dk,i_State);
        }
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
                                     Conservative_s[c],
                                     Conservative_s[xm],
#if (AMREX_SPACEDIM > 1)
                                     Conservative_s[xp], Conservative_s[ym], Conservative_s[yp],
                                     Conservative_s[xm_ym], Conservative_s[xm_yp], Conservative_s[xp_ym],
#endif
#if (AMREX_SPACEDIM > 2)
                                     Conservative_s[zm], Conservative_s[zp],
                                     Conservative_s[xm_zm], Conservative_s[xm_zp], Conservative_s[xp_zm],
                                     Conservative_s[ym_zm], Conservative_s[ym_zp], Conservative_s[yp_zm],
#endif
                                     mu, mub, Kappa,
                        AMREX_D_DECL(dx, dy, dz),
                                     EOS_Parameter,
                                     ID_Dimension);
    
}

#endif
#endif