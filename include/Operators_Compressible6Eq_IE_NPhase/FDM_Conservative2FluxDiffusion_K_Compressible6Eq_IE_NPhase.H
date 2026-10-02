// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_K_Compressible6Eq_IE_NPhase_H_
#define FDM_Conservative2FluxDiffusion_K_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// Viscous and conductive flux on the lower face of cell (i,j,k) in direction ID_Dimension.
// Gathers the face stencil (FaceStencil::OffsetOnFace, 1 ghost cell with the edge
// neighbors) and calls FDM_Conservative2FluxDiffusion_D.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_K(int i, int j, int k,
                                      amrex::Array<amrex::Real,NSTATE>& Flux,
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
    amrex::ignore_unused(prob_lo, time);

    amrex::Array<amrex::Array<amrex::Real,NSTATE>,FaceStencil::N> Conservative_s;
    for (int s = 0; s < FaceStencil::N; ++s){
        int di, dj, dk;
        FaceStencil::OffsetOnFace(s, ID_Dimension, di, dj, dk);
        for (int iState = 0; iState < NSTATE; ++iState){
            Conservative_s[s][iState] = Conservative(i+di, j+dj, k+dk, iState);
        }
    }

    FDM_Conservative2FluxDiffusion_D(Flux,
                                     dmax,
                                     Conservative_s,
                        AMREX_D_DECL(dx, dy, dz),
                                     EOS_Parameter,
                                     ID_Dimension);
}

#endif
#endif
