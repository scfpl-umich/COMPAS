// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_D_Compressible6Eq_H_
#define FDM_Conservative2FluxDiffusion_D_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS && DIFFUSION == true)

// =======================================
// From the conservative states of the face stencil (Conservative_s[s], s a FaceStencil point)
// to the face quantities of FluxDiffusion:
//   - the face velocity and velocity gradient (FDM_FaceVelocityGradient_D);
//   - the face viscosities, the harmonic mean of the mixture viscosities of the two cells
//     (PhaseDiffusion_FaceViscosity_D);
//   - per phase, the temperature difference across the face, the harmonic face weight and
//     the thermal diffusivity bound for the time step (where the phase conducts);
//   - with the pressure-temperature relaxation, the heat flux of the mixture instead (below).
//
// Pressure-temperature relaxation. After the relaxation the phases of a cell have one
// temperature, and the relaxed state depends only on the phase masses, the momentum and the
// total energy, not on how the energy is shared between the phases (Huang, arXiv:2609.18085,
// Sec. 3.1.2). The phases then conduct as one medium, q = -kappa grad T with the mixture
// conductivity kappa = alpha_1 kappa_1 + alpha_2 kappa_2 (the sum of the phase heat fluxes),
// discretized with the harmonic mean of the conductivities of the two cells (as the viscosity,
// PhaseDiffusion_FaceViscosity_D) and the cell temperature, the mean of the temperatures of the
// phases present weighted by their heat capacities alpha_k rho_k cv_k: the temperature the cell
// reaches when its phases exchange heat at fixed volumes (de_k = cv_k dT_k for every EOS with a
// temperature). After the relaxation the phase temperatures are equal and any weighting gives
// that temperature; the weighting matters where the state is not relaxed (a cell where the
// pressure-temperature relaxation fell back to the pressure relaxation, interpolated coarse-fine
// ghost cells, a coarse cell after the reflux). In a uniform mixture this is the sum of the phase heat
// fluxes; across a material interface it is the series conductance of the two cells, which
// the phase heat fluxes cannot carry (a phase conducts only through itself, so a sharp interface
// would not conduct). The flux goes into Q_Mix, which FluxDiffusion adds to the energy of phase
// 1 and to the work component, so FVM_SurfaceIntegral_NC gives each phase its mass fraction
// of the heat. The diffusivity bound uses the heat capacity of the mixture, the smaller
// sum_k alpha_k rho_k cv_k of the two cells.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_D(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& Flux,
                                      amrex::Real& Flux_Work,
                                      amrex::Real& dmax,
                                      amrex::Array<amrex::Array<amrex::Real,5+AMREX_SPACEDIM>,FaceStencil::N> const& Conservative_s,
                         AMREX_D_DECL(amrex::Real const& dx,
                                      amrex::Real const& dy,
                                      amrex::Real const& dz),
                                      Parm const& EOS_Parameter,
                                      int ID_Dimension)
{
    // Velocity at the stencil cells, then on the face
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
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Grad_Vel;
    FDM_FaceVelocityGradient_D(Vel,
                               Grad_Vel,
                               V,
                  AMREX_D_DECL(dx, dy, dz),
                               ID_Dimension);

    // The two cells across the face
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative_m = Conservative_s[FaceStencil::Minus(ID_Dimension)];
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative_c = Conservative_s[FaceStencil::c];

    // Face viscosities
    amrex::Real mu_m, mub_m, mu_c, mub_c;
    Conservative2Viscosity(mu_m, mub_m, Conservative_m, EOS_Parameter);
    Conservative2Viscosity(mu_c, mub_c, Conservative_c, EOS_Parameter);
    amrex::Real mu, mub;
    PhaseDiffusion_FaceViscosity_D(mu, mu_m, mu_c);
    PhaseDiffusion_FaceViscosity_D(mub, mub_m, mub_c);

    // Phase temperatures, face weights and temperature gradients
    amrex::Array<amrex::Real,2> T_m, T_c;
    amrex::Array<bool,2> Present_m, Present_c;
    Conservative2PhaseTemperature(T_m, Present_m, Conservative_m, EOS_Parameter);
    Conservative2PhaseTemperature(T_c, Present_c, Conservative_c, EOS_Parameter);

    amrex::Array<amrex::Real,2> const Alpha_m     = {Conservative_m[INDEX_VolumeFraction1], 1.0 - Conservative_m[INDEX_VolumeFraction1]};
    amrex::Array<amrex::Real,2> const Alpha_c     = {Conservative_c[INDEX_VolumeFraction1], 1.0 - Conservative_c[INDEX_VolumeFraction1]};
    amrex::Array<amrex::Real,2> const Alpha_rho_m = {Conservative_m[INDEX_Mass1], Conservative_m[INDEX_Mass2]};
    amrex::Array<amrex::Real,2> const Alpha_rho_c = {Conservative_c[INDEX_Mass1], Conservative_c[INDEX_Mass2]};
    amrex::Array<amrex::Real,2> const Kappa       = {EOS_Parameter.kappa_1, EOS_Parameter.kappa_2};
    amrex::Array<amrex::Real,2> const Cv          = {EOS_Parameter.cv_1, EOS_Parameter.cv_2};
    amrex::Array<amrex::Real,AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx, dy, dz)};

    amrex::Array<amrex::Real,2> Weight = {0.0, 0.0};
    amrex::Array<amrex::Real,2> Grad_T = {0.0, 0.0};
    amrex::Real Chi_T = 0.0;
    amrex::Real Q_Mix = 0.0;
    if (EOS_Parameter.Physics_Parm.pressure_temperature_relaxation){
        // The mixture conducts as one medium (see above), in place of the phase heat fluxes
        amrex::Real Kappa_m = Alpha_m[0]*Kappa[0] + Alpha_m[1]*Kappa[1];
        amrex::Real Kappa_c = Alpha_c[0]*Kappa[0] + Alpha_c[1]*Kappa[1];
        amrex::Real Kappa_Face;
        PhaseDiffusion_FaceViscosity_D(Kappa_Face, Kappa_m, Kappa_c);
        amrex::Real rhoCv_m = 0.0, rhoCv_c = 0.0, rhoCvT_m = 0.0, rhoCvT_c = 0.0;
        for (int kp = 0; kp < 2; ++kp){
            if (Present_m[kp]){
                rhoCv_m  += Alpha_rho_m[kp]*Cv[kp];
                rhoCvT_m += Alpha_rho_m[kp]*Cv[kp]*T_m[kp];
            }
            if (Present_c[kp]){
                rhoCv_c  += Alpha_rho_c[kp]*Cv[kp];
                rhoCvT_c += Alpha_rho_c[kp]*Cv[kp]*T_c[kp];
            }
        }
        if (Kappa_Face > 0.0 && rhoCv_m > 0.0 && rhoCv_c > 0.0){
            amrex::Real Grad_T_Mix;
            FDM_Gradient_D(Grad_T_Mix, rhoCvT_m/rhoCv_m, rhoCvT_c/rhoCv_c, dX[ID_Dimension]);
            Q_Mix = Kappa_Face*Grad_T_Mix;
            Chi_T = Kappa_Face/std::min(rhoCv_m, rhoCv_c);
        }
    }
    else {
        // Each phase conducts through itself
        for (int kp = 0; kp < 2; ++kp){
            PhaseDiffusion_FaceWeight_D(Weight[kp], Alpha_m[kp], Alpha_c[kp], Present_m[kp], Present_c[kp]);
            if (Weight[kp] > 0.0){
                FDM_Gradient_D(Grad_T[kp], T_m[kp], T_c[kp], dX[ID_Dimension]);
                // The face weight over the smaller volume fraction (between 1 and 2) bounds the
                // explicit limit of both cells
                amrex::Real rho_k = std::min(Alpha_rho_m[kp]/Alpha_m[kp], Alpha_rho_c[kp]/Alpha_c[kp]);
                amrex::Real Ratio = Weight[kp]/std::min(Alpha_m[kp], Alpha_c[kp]);
                Chi_T = std::max(Chi_T, Ratio*Kappa[kp]/(rho_k*Cv[kp]));
            }
        }
    }

    amrex::Real rho_m = Alpha_rho_m[0] + Alpha_rho_m[1];
    amrex::Real rho_c = Alpha_rho_c[0] + Alpha_rho_c[1];
    amrex::Real rho   = std::min(rho_m, rho_c);

    FluxDiffusion(Flux,
                  Flux_Work,
                  dmax,
                  Vel,
                  Grad_Vel,
                  Weight,
                  Grad_T,
                  Q_Mix,
                  mu,
                  mub,
                  rho,
                  Chi_T,
                  EOS_Parameter,
                  ID_Dimension);
}

#endif
#endif
