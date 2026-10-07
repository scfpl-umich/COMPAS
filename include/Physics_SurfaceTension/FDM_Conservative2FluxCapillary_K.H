// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxCapillary_K_H_
#define FDM_Conservative2FluxCapillary_K_H_

#include<Macros.H>
#if (SURFACE_TENSION == true)

// =======================================
// Capillary stress on the lower face of cell (i,j,k) in direction ID_Dimension, for every model.
//
//   T.e_d = sum_k sigma_k ( |g_k| e_d - g_kd g_k/|g_k| )      (0 where |g_k| is tiny)
//
// g_k is the compact face gradient of alpha_k (normal derivative from the two cells across the
// face, tangential ones from the edge averages of four cells: 1 ghost cell with corners), from
// FDM_Gradient_{X,Y,Z}FaceFromCenter_K. The two-phase models use T.e_d = sigma_12 (|g| e_d - g_d g/|g|)
// with g from alpha_1. The N-phase models take g_k of the stored alpha_k for the first N - 1 phases
// and g_N = -sum of the others, so that sum_k g_k = 0 on every face and N = 2 reduces to the
// two-phase stress. For a uniform interface profile, div T = sum_k sigma_k kappa_k grad alpha_k
// with kappa_k = -div(g_k/|g_k|).
//
// The volume fractions come from a separate array A (component n holds alpha_{n+1}, ST_NALPHA
// components), which holds the stored ones smoothed SurfaceTension.smoothing times by the 1-2-1
// filter (SurfaceTension_Smooth_K; a copy without smoothing).
//
// The flux follows the sign of the diffusive flux (U_t = -div(F - F_ST)):
//   F_ST = [0, ..., T.e_d (momentum), u_f.T.e_d (energy, energy_form = conservative only), ..., 0]
// with u_f = (u_m + u_c)/2. The energy is the mixture energy of the five-equation models, rho E of
// SIXEQS_IE_NPHASE, and energy 1 of SIXEQS (FVM_SurfaceIntegral_NC_Capillary moves the share of
// phase 2). The NC components INDEX_NC_CapillaryStress + b hold T.e_d, and INDEX_NC_CapillaryWork
// (SIXEQS) holds u_f.T.e_d.
// =======================================

#if (PHYSICS == SIXEQS)
inline constexpr int INDEX_Capillary_Energy = INDEX_Energy1;
#elif (PHYSICS == SIXEQS_IE_NPHASE)
inline constexpr int INDEX_Capillary_Energy = INDEX_Total_Energy;
#else
inline constexpr int INDEX_Capillary_Energy = INDEX_Energy;
#endif

// Below this |g| h (h the smallest cell size) the stress of a phase is 0
inline constexpr amrex::Real Capillary_Gradient_Tiny = 1.0e-12;

// Mixture density of cell (i,j,k)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real SurfaceTension_Density(int i, int j, int k,
                                   amrex::Array4<amrex::Real const> const& U)
{
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    return U(i,j,k,INDEX_Mass1) + U(i,j,k,INDEX_Mass2);
#else
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += U(i,j,k,INDEX_Mass1+i_Phase);
    }
    return rho;
#endif
}

// Volume fraction of phase i_Phase (0-based) of cell (i,j,k)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real SurfaceTension_VolumeFraction(int i, int j, int k, int i_Phase,
                                          amrex::Array4<amrex::Real const> const& U)
{
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    return (i_Phase == 0) ? U(i,j,k,INDEX_VolumeFraction1) : 1.0 - U(i,j,k,INDEX_VolumeFraction1);
#else
    return U(i,j,k,INDEX_VolumeFraction1+i_Phase);
#endif
}

// Velocity of cell (i,j,k)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SurfaceTension_Velocity(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Vel,
                             int i, int j, int k,
                             amrex::Array4<amrex::Real const> const& U)
{
    amrex::Real rho = SurfaceTension_Density(i, j, k, U);
    for (int b = 0; b < AMREX_SPACEDIM; b++){
        Vel[b] = U(i,j,k,INDEX_MomentumX+b)/rho;
    }
}

// Face gradient of component n on the lower face of (i,j,k) in direction ID_Dimension
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SurfaceTension_FaceGradient_K(int i, int j, int k, int n,
                                   amrex::Array<amrex::Real,AMREX_SPACEDIM>& Grad,
                                   amrex::Array4<amrex::Real const> const& U,
                      AMREX_D_DECL(amrex::Real dx,
                                   amrex::Real dy,
                                   amrex::Real dz),
                                   int ID_Dimension)
{
    if (ID_Dimension == 0){
        FDM_Gradient_XFaceFromCenter_K(i, j, k, n, Grad, U, AMREX_D_DECL(dx, dy, dz));
    }
#if (AMREX_SPACEDIM > 1)
    else if (ID_Dimension == 1){
        FDM_Gradient_YFaceFromCenter_K(i, j, k, n, Grad, U, AMREX_D_DECL(dx, dy, dz));
    }
#endif
#if (AMREX_SPACEDIM > 2)
    else if (ID_Dimension == 2){
        FDM_Gradient_ZFaceFromCenter_K(i, j, k, n, Grad, U, AMREX_D_DECL(dx, dy, dz));
    }
#endif
}

// One pass of the 1-2-1 filter (the tensor product of the weights 1/4, 1/2, 1/4 in each direction)
// on component n of In; reads one layer of neighbors, corners included
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real SurfaceTension_Smooth_K(int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& In)
{
    constexpr amrex::Real w[3] = {0.25, 0.5, 0.25};
    amrex::Real Sum = 0.0;
#if (AMREX_SPACEDIM > 2)
    for (int c = -1; c <= 1; c++){
#else
    {
        const int c = 0;
#endif
#if (AMREX_SPACEDIM > 1)
        for (int b = -1; b <= 1; b++){
#else
        {
            const int b = 0;
#endif
            for (int a = -1; a <= 1; a++){
                amrex::Real W = w[a+1];
#if (AMREX_SPACEDIM > 1)
                W *= w[b+1];
#endif
#if (AMREX_SPACEDIM > 2)
                W *= w[c+1];
#endif
                Sum += W*In(i+a, j+b, k+c, n);
            }
        }
    }
    return Sum;
}

// Stress += Sigma ( |g| e_d - g_d g/|g| )
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void CapillaryStress_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Stress,
                       amrex::Array<amrex::Real,AMREX_SPACEDIM> const& Grad,
                       amrex::Real Sigma,
                       amrex::Real Tiny,
                       int ID_Dimension)
{
    amrex::Real Norm = std::sqrt(AMREX_D_TERM(Grad[0]*Grad[0], + Grad[1]*Grad[1], + Grad[2]*Grad[2]));
    if (Norm > Tiny){
        for (int b = 0; b < AMREX_SPACEDIM; b++){
            amrex::Real Delta = (b == ID_Dimension) ? Norm : 0.0;
            Stress[b] += Sigma*(Delta - Grad[ID_Dimension]*Grad[b]/Norm);
        }
    }
}

// Capillary stress T.e_d on the lower face of (i,j,k) in direction ID_Dimension, from the volume
// fractions A (component n: alpha_{n+1})
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_CapillaryStress_K(int i, int j, int k,
                           amrex::Array<amrex::Real,AMREX_SPACEDIM>& Stress,
                           amrex::Array4<amrex::Real const> const& A,
              AMREX_D_DECL(amrex::Real dx,
                           amrex::Real dy,
                           amrex::Real dz),
                           int ID_Dimension,
                           SurfaceTension_Parameter const& ST)
{
    for (int b = 0; b < AMREX_SPACEDIM; b++){
        Stress[b] = 0.0;
    }
#if (AMREX_SPACEDIM == 1)
    amrex::Real h = dx;
#else
    amrex::Real h = amrex::min(AMREX_D_DECL(dx, dy, dz));
#endif
    amrex::Real Tiny = Capillary_Gradient_Tiny/h;

    amrex::Array<amrex::Real,AMREX_SPACEDIM> Grad;
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    SurfaceTension_FaceGradient_K(i, j, k, 0, Grad, A, AMREX_D_DECL(dx, dy, dz), ID_Dimension);
    CapillaryStress_D(Stress, Grad, ST.sigma_12, Tiny, ID_Dimension);
#else
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Grad_Last;
    for (int b = 0; b < AMREX_SPACEDIM; b++){
        Grad_Last[b] = 0.0;
    }
    for (int i_Phase = 0; i_Phase < NPHASE-1; i_Phase++){
        SurfaceTension_FaceGradient_K(i, j, k, i_Phase, Grad, A, AMREX_D_DECL(dx, dy, dz), ID_Dimension);
        CapillaryStress_D(Stress, Grad, ST.sigma[i_Phase], Tiny, ID_Dimension);
        for (int b = 0; b < AMREX_SPACEDIM; b++){
            Grad_Last[b] -= Grad[b];
        }
    }
    CapillaryStress_D(Stress, Grad_Last, ST.sigma[NPHASE-1], Tiny, ID_Dimension);
#endif
}

// Capillary face flux F_ST and the NC face components
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxCapillary_K(int i, int j, int k,
                                      amrex::Array<amrex::Real,NSTATE>& Flux,
                                      amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
                                      amrex::Array4<amrex::Real const> const& U,
                                      amrex::Array4<amrex::Real const> const& A,
                         AMREX_D_DECL(amrex::Real dx,
                                      amrex::Real dy,
                                      amrex::Real dz),
                                      int ID_Dimension,
                                      SurfaceTension_Parameter const& ST)
{
    for (int iState = 0; iState < NSTATE; iState++){
        Flux[iState] = 0.0;
    }

    amrex::Array<amrex::Real,AMREX_SPACEDIM> Stress;
    FDM_CapillaryStress_K(i, j, k, Stress, A, AMREX_D_DECL(dx, dy, dz), ID_Dimension, ST);

    // Face velocity u_f = (u_m + u_c)/2
    int im = (ID_Dimension == 0) ? i-1 : i;
    int jm = (ID_Dimension == 1) ? j-1 : j;
    int km = (ID_Dimension == 2) ? k-1 : k;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel_m, Vel_c;
    SurfaceTension_Velocity(Vel_m, im, jm, km, U);
    SurfaceTension_Velocity(Vel_c, i, j, k, U);

    amrex::Real Work = 0.0;
    for (int b = 0; b < AMREX_SPACEDIM; b++){
        amrex::Real Vel_f;
        Interpolation2Face_Linear_D(Vel_f, Vel_m[b], Vel_c[b]);
        Work += Vel_f*Stress[b];
        Flux[INDEX_MomentumX+b] = Stress[b];
        FluxNC[INDEX_NC_CapillaryStress+b] = Stress[b];
    }
    if (ST.energy_form == CapillaryEnergy::Conservative){
        Flux[INDEX_Capillary_Energy] = Work;
    }
#if (PHYSICS == SIXEQS)
    FluxNC[INDEX_NC_CapillaryWork] = Work;
#endif
}

// =======================================
// Capillary energy of a cell from the divergences of the NC components (dt-weighted increments in
// the reflux): Div_Stress = [div T]_h, Div_Work = [div(u.T)]_h (SIXEQS), with the velocity and
// mass fractions of the state U at (i,j,k).
//   energy_form = work: the energy gets u_c.[div T]_h; SIXEQS splits it as Y_k u_c.[div T]_h.
//   energy_form = conservative: the energy flux already holds u_f.T.e_d; SIXEQS moves
//   Y_2 [div(u.T)]_h from energy 1 (which received all of it) to energy 2.
// =======================================
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SurfaceTension_Energy_D(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& dof_out,
                             amrex::Array4<amrex::Real const> const& U,
                             amrex::Array<amrex::Real,AMREX_SPACEDIM> const& Div_Stress,
                             amrex::Real Div_Work,
                             SurfaceTension_Parameter const& ST)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    SurfaceTension_Velocity(Vel, i, j, k, U);
#if (PHYSICS == SIXEQS)
    amrex::Real Y2 = U(i,j,k,INDEX_Mass2)/SurfaceTension_Density(i, j, k, U);
#endif

    if (ST.energy_form == CapillaryEnergy::Work){
        amrex::Real Power = 0.0;
        for (int b = 0; b < AMREX_SPACEDIM; b++){
            Power += Vel[b]*Div_Stress[b];
        }
        dof_out(i,j,k,INDEX_Capillary_Energy) += Power;
#if (PHYSICS == SIXEQS)
        amrex::Real Transfer = Y2*Power;
        dof_out(i,j,k,INDEX_Energy1) -= Transfer;
        dof_out(i,j,k,INDEX_Energy2) += Transfer;
#endif
    }
#if (PHYSICS == SIXEQS)
    else {
        amrex::Real Transfer = Y2*Div_Work;
        dof_out(i,j,k,INDEX_Energy1) -= Transfer;
        dof_out(i,j,k,INDEX_Energy2) += Transfer;
    }
#else
    amrex::ignore_unused(Div_Work);
#endif
}

// Capillary NC terms of the cell update (after FVM_SurfaceIntegral_NC)
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_SurfaceIntegral_NC_Capillary(int i, int j, int k,
                                      amrex::Array4<amrex::Real> const& dof_out,
                                      amrex::Array4<amrex::Real const> const& dof_in,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxxNC,
                                      amrex::Array4<amrex::Real const> const& flxyNC,
                                      amrex::Array4<amrex::Real const> const& flxzNC),
                         AMREX_D_DECL(amrex::Real dx,
                                      amrex::Real dy,
                                      amrex::Real dz),
                                      SurfaceTension_Parameter const& ST)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Div_Stress;
    for (int b = 0; b < AMREX_SPACEDIM; b++){
        FVM_Divergence_K(i, j, k, INDEX_NC_CapillaryStress+b,
                         Div_Stress[b],
            AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
            AMREX_D_DECL(dx    , dy    , dz    ));
    }
    amrex::Real Div_Work = 0.0;
#if (PHYSICS == SIXEQS)
    FVM_Divergence_K(i, j, k, INDEX_NC_CapillaryWork,
                     Div_Work,
        AMREX_D_DECL(flxxNC, flxyNC, flxzNC),
        AMREX_D_DECL(dx    , dy    , dz    ));
#endif
    SurfaceTension_Energy_D(i, j, k, dof_out, dof_in, Div_Stress, Div_Work, ST);
}

#endif
#endif
