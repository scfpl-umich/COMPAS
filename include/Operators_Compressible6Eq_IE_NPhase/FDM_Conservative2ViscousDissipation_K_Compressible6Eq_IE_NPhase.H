// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2ViscousDissipation_K_Compressible6Eq_IE_NPhase_H_
#define FDM_Conservative2ViscousDissipation_K_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// Viscous dissipation Phi_k = tau_k : grad u of each phase at the center of cell (i,j,k)
// (ViscousDissipation_D), with grad u from central differences of u = rho u/rho over the
// face neighbors (1 ghost cell).
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2ViscousDissipation_K(int i, int j, int k,
                                           amrex::Array<amrex::Real,NPHASE>& Phi,
                                           amrex::Array4<amrex::Real const> const& Conservative,
                              AMREX_D_DECL(amrex::Real const& dx,
                                           amrex::Real const& dy,
                                           amrex::Real const& dz),
                                           Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Vel_m, Vel_p;
    for (int a = 0; a < AMREX_SPACEDIM; ++a){
        int const di = (a == 0), dj = (a == 1), dk = (a == 2);
        amrex::Array<amrex::Real,NSTATE> Conservative_m, Conservative_p;
        for (int iState = 0; iState < NSTATE; ++iState){
            Conservative_m[iState] = Conservative(i-di, j-dj, k-dk, iState);
            Conservative_p[iState] = Conservative(i+di, j+dj, k+dk, iState);
        }
        AMREX_D_TERM( amrex::Real u;,
                      amrex::Real v;,
                      amrex::Real w );
        Conservative2Velocity(AMREX_D_DECL(u, v, w), Conservative_m);
        AMREX_D_TERM( Vel_m[a][0] = u;,
                      Vel_m[a][1] = v;,
                      Vel_m[a][2] = w );
        Conservative2Velocity(AMREX_D_DECL(u, v, w), Conservative_p);
        AMREX_D_TERM( Vel_p[a][0] = u;,
                      Vel_p[a][1] = v;,
                      Vel_p[a][2] = w );
    }

    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Grad_Vel;
    FDM_CellVelocityGradient_D(Grad_Vel,
                               Vel_m,
                               Vel_p,
                  AMREX_D_DECL(dx, dy, dz));

    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        ViscousDissipation_D(Phi[i_Phase], Grad_Vel, EOS_Parameter.mu[i_Phase], EOS_Parameter.muB[i_Phase]);
    }
}

#endif
#endif
