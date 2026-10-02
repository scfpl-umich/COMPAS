// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_
#define FDM_Conservative2FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// From the conservative states of the face stencil (Conservative_s[s], s a FaceStencil point)
// to the cell values that FDM_FluxDiffusion_D needs: the velocity of every stencil cell, and
// for the two cells across the face the mixture viscosities, the density and, per phase, the
// volume fraction, the partial density and the temperature where the phase is present.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_D(amrex::Array<amrex::Real,NSTATE>& Flux,
                                      amrex::Real& dmax,
                                      amrex::Array<amrex::Array<amrex::Real,NSTATE>,FaceStencil::N> const& Conservative_s,
                         AMREX_D_DECL(amrex::Real const& dx,
                                      amrex::Real const& dy,
                                      amrex::Real const& dz),
                                      Parm const& EOS_Parameter,
                                      int ID_Dimension)
{
    // Velocity at the stencil cells
    amrex::Array<amrex::Array<amrex::Real,FaceStencil::N>,AMREX_SPACEDIM> V;
    for (int s = 0; s < FaceStencil::N; ++s){
        AMREX_D_TERM( amrex::Real u;,
                      amrex::Real v;,
                      amrex::Real w );
        Conservative2Velocity(AMREX_D_DECL(u, v, w),
                              Conservative_s[s]);
        AMREX_D_TERM( V[0][s] = u;,
                      V[1][s] = v;,
                      V[2][s] = w );
    }

    // The two cells across the face
    amrex::Array<amrex::Real,NSTATE> const& Conservative_m = Conservative_s[FaceStencil::Minus(ID_Dimension)];
    amrex::Array<amrex::Real,NSTATE> const& Conservative_c = Conservative_s[FaceStencil::c];

    amrex::Real mu_m, mub_m, mu_c, mub_c;
    Conservative2Viscosity(mu_m, mub_m, Conservative_m, EOS_Parameter);
    Conservative2Viscosity(mu_c, mub_c, Conservative_c, EOS_Parameter);

    amrex::Array<amrex::Real,NPHASE> T_m, T_c;
    amrex::Array<bool,NPHASE> Present_m, Present_c;
    Conservative2PhaseTemperature(T_m, Present_m, Conservative_m, EOS_Parameter);
    Conservative2PhaseTemperature(T_c, Present_c, Conservative_c, EOS_Parameter);

    amrex::Array<amrex::Real,NPHASE> Alpha_m, Alpha_c, Alpha_rho_m, Alpha_rho_c;
    amrex::Real rho_m = 0.0;
    amrex::Real rho_c = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        Alpha_m[i_Phase]     = Conservative_m[INDEX_VolumeFraction1 + i_Phase];
        Alpha_c[i_Phase]     = Conservative_c[INDEX_VolumeFraction1 + i_Phase];
        Alpha_rho_m[i_Phase] = Conservative_m[INDEX_Mass1 + i_Phase];
        Alpha_rho_c[i_Phase] = Conservative_c[INDEX_Mass1 + i_Phase];
        rho_m += Alpha_rho_m[i_Phase];
        rho_c += Alpha_rho_c[i_Phase];
    }

    FDM_FluxDiffusion_D(Flux,
                        dmax,
                        V,
                        mu_m, mu_c,
                        mub_m, mub_c,
                        rho_m, rho_c,
                        Alpha_m, Alpha_c,
                        Alpha_rho_m, Alpha_rho_c,
                        T_m, T_c,
                        Present_m, Present_c,
           AMREX_D_DECL(dx, dy, dz),
                        EOS_Parameter,
                        ID_Dimension);
}

#endif
#endif
