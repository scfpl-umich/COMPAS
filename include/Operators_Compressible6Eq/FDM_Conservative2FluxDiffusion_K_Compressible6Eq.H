// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_K_Compressible6Eq_H_
#define FDM_Conservative2FluxDiffusion_K_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS && DIFFUSION == true)

// =======================================
// Viscous and conductive flux on the lower face of cell (i,j,k) in direction ID_Dimension.
// Gathers the face stencil (FaceStencil::OffsetOnFace, 1 ghost cell with the edge
// neighbors) and calls FDM_Conservative2FluxDiffusion_D. Besides the conservative flux, it
// writes the viscous energy flux (u.tau)_n (plus the heat flux of the mixture with the
// pressure-temperature relaxation) to the non-conservative component
// INDEX_NC_ViscousWork; the other components of FluxNC are left as the hyperbolic kernel set
// them, so it runs after it.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_K(int i, int j, int k,
                                      amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Flux,
                                      amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
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

    amrex::Array<amrex::Array<amrex::Real,5+AMREX_SPACEDIM>,FaceStencil::N> Conservative_s;
    for (int s = 0; s < FaceStencil::N; ++s){
        int di, dj, dk;
        FaceStencil::OffsetOnFace(s, ID_Dimension, di, dj, dk);
        for (int iState = 0; iState < 5+AMREX_SPACEDIM; ++iState){
            Conservative_s[s][iState] = Conservative(i+di, j+dj, k+dk, iState);
        }
    }

    FDM_Conservative2FluxDiffusion_D(Flux,
                                     FluxNC[INDEX_NC_ViscousWork],
                                     dmax,
                                     Conservative_s,
                        AMREX_D_DECL(dx, dy, dz),
                                     EOS_Parameter,
                                     ID_Dimension);
}

#endif
#endif
