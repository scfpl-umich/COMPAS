// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_WENO5_Quad_K_Compressible5Eq_H_
#define FVM_Reconstruction2FaceLR_Conservative_WENO5_Quad_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

#include <cmath>
#include <array>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void WENO5_Quad(
    amrex::Array<amrex::Real,3>& f_WENO,
    amrex::Array<amrex::Real,3>& Weight_WENO,
    amrex::Real const& f_mm,
    amrex::Real const& f_m,
    amrex::Real const& f_c,
    amrex::Real const& f_p,
    amrex::Real const& f_pp,
    int ID_Node)
{
    constexpr amrex::Real sq3 = 1.7320508075688772;

    amrex::Real d0, d1, d2;

    // =====================================================
    // Candidate polynomials
    // =====================================================

    // --------------------
    // Right face x_{i+1/2}
    // --------------------
    if (ID_Node == 1)
    {
        f_WENO[0] = f_c + (2.0/3.0)*(f_p - f_c) - (1.0/6.0)*(f_pp - f_p);
        f_WENO[1] = f_c + (1.0/6.0)*(f_c - f_m) + (1.0/3.0)*(f_p - f_c);
        f_WENO[2] = f_c - (1.0/3.0)*(f_m - f_mm) + (5.0/6.0)*(f_c - f_m);
              
        d0 = 3.0/10.0;
        d2 = 1.0/10.0;
        d1 = 1.0 - d0 - d2;
    }

    // --------------------
    // Left face x_{i-1/2}
    // --------------------
    else if (ID_Node == -1)
    {
        f_WENO[0] = f_c - (5.0/6.0)*(f_p - f_c) + (1.0/3.0)*(f_pp - f_p);
        f_WENO[1] = f_c - (1.0/3.0)*(f_c - f_m) - (1.0/6.0)*(f_p - f_c);
        f_WENO[2] = f_c + (1.0/6.0)*(f_m - f_mm) - (2.0/3.0)*(f_c - f_m);
              
        d0 = 1.0/10.0;
        d2 = 3.0/10.0;
        d1 = 1.0 - d0 - d2;
    }

    // --------------------
    // Gauss left x_i - dx/(2√3)
    // --------------------
    else if (ID_Node == 2)
    {
        f_WENO[0] = f_c + (3*f_c - 4*f_p + f_pp)*( sq3 / 12.0 );
        f_WENO[1] = f_c - (-f_m + f_p)*( sq3 / 12.0 );
        f_WENO[2] = f_c - (3*f_c - 4*f_m + f_mm)*( sq3 / 12.0 );

        d0 = (210.0 - sq3)/1080.0;
        d1 = 11.0/18.0;
        d2 = (210.0 + sq3)/1080.0;
    }

    // --------------------
    // Gauss right x_i + dx/(2√3)
    // --------------------
    else if (ID_Node == 3)
    {
        f_WENO[0] = f_c - (3*f_c - 4*f_p + f_pp)*( sq3 / 12.0 );
        f_WENO[1] = f_c - (f_m - f_p)*( sq3 / 12.0 );
        f_WENO[2] = f_c - (-3*f_c + 4*f_m - f_mm)*( sq3 / 12.0 );

        d0 = (210.0 + sq3)/1080.0;
        d1 = 11.0/18.0;
        d2 = (210.0 - sq3)/1080.0;
    }

    // =====================================================
    // Smoothness indicators (same for all nodes)
    // =====================================================

    amrex::Real IS0 =
        (13./12.)*(f_c - 2*f_p + f_pp)*(f_c - 2*f_p + f_pp)
      + 0.25*(3*f_c - 4*f_p + f_pp)*(3*f_c - 4*f_p + f_pp);

    amrex::Real IS1 =
        (13./12.)*(f_m - 2*f_c + f_p)*(f_m - 2*f_c + f_p)
      + 0.25*(f_m - f_p)*(f_m - f_p);

    amrex::Real IS2 =
        (13./12.)*(f_mm - 2*f_m + f_c)*(f_mm - 2*f_m + f_c)
      + 0.25*(f_mm - 4*f_m + 3*f_c)*(f_mm - 4*f_m + 3*f_c);

    // =====================================================
    // Ratio-based nonlinear weights (no epsilon)
    // =====================================================

    // amrex::Real r01 = IS0/IS1; if (std::isnan(r01)) r01 = 1.0;
    // amrex::Real r10 = 1.0/r01;

    // amrex::Real r02 = IS0/IS2; if (std::isnan(r02)) r02 = 1.0;
    // amrex::Real r20 = 1.0/r02;

    // amrex::Real r12 = IS1/IS2; if (std::isnan(r12)) r12 = 1.0;
    // amrex::Real r21 = 1.0/r12;

    // Weight_WENO[0] =
    //     d0 / ( d0 + (r01*r01)*d1 + (r02*r02)*d2 );

    // Weight_WENO[1] =
    //     d1 / ( (r10*r10)*d0 + d1 + (r12*r12)*d2 );

    // Weight_WENO[2] =
    //     d2 / ( (r20*r20)*d0 + (r21*r21)*d1 + d2 );



    // amrex::Real eps = 1.0e-40; 
    // amrex::Real a0 = d0 / ( (IS0 + eps)*(IS0 + eps) );
    // amrex::Real a1 = d1 / ( (IS1 + eps)*(IS1 + eps) );
    // amrex::Real a2 = d2 / ( (IS2 + eps)*(IS2 + eps) );

    // amrex::Real sum_a = a0 + a1 + a2;

    // Weight_WENO[0] = a0 / sum_a;
    // Weight_WENO[1] = a1 / sum_a;
    // Weight_WENO[2] = a2 / sum_a;




    // Standard WENO-Z with proper epsilon
    amrex::Real tau5 = std::abs(IS0 - IS2);
    amrex::Real eps = 1.0e-6; // Buries double-precision noise

    // Power p=1 is sufficient and less prone to oscillatory noise
    amrex::Real a0 = d0 * (1.0 + (tau5 / (IS0 + eps)));
    amrex::Real a1 = d1 * (1.0 + (tau5 / (IS1 + eps)));
    amrex::Real a2 = d2 * (1.0 + (tau5 / (IS2 + eps)));

    amrex::Real sum_a = a0 + a1 + a2;

    Weight_WENO[0] = a0 / sum_a;
    Weight_WENO[1] = a1 / sum_a;
    Weight_WENO[2] = a2 / sum_a;
}
// 1e-1
// Final Convergence Rates (last two points):
// -------------------------------------------------
// a1rho1               :   4.6892
// a2rho2               : ERROR computing rate (Non-positive error encountered)
// rho_u                :   4.3730
// rho_v                :   4.2474
// rho_E                :   4.8735
// a1                   :  -0.0000
// -------------------------------------------------

// 1e-3
// Final Convergence Rates (last two points):
// -------------------------------------------------
// a1rho1               :   4.6913
// a2rho2               : ERROR computing rate (Non-positive error encountered)
// rho_u                :   4.3798
// rho_v                :   4.2539
// rho_E                :   4.8750
// a1                   :  -0.0000
// -------------------------------------------------

// 1e-6
// Final Convergence Rates (last two points):
// -------------------------------------------------
// a1rho1               :   4.6099
// a2rho2               : ERROR computing rate (Non-positive error encountered)
// rho_u                :   4.3528
// rho_v                :   4.2075
// rho_E                :   4.8922
// a1                   :  -0.0000
// -------------------------------------------------


// 1e-12
// Final Convergence Rates (last two points):
// -------------------------------------------------
// a1rho1               :   4.1696
// a2rho2               : ERROR computing rate (Non-positive error encountered)
// rho_u                :   4.4000
// rho_v                :   4.2116
// rho_E                :   4.6229
// a1                   :  -0.0000
// -------------------------------------------------

template<typename T_WENO>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Reconstruction2Node_WENO(
    amrex::Real& f_Node,
    T_WENO const& f_WENO,
    T_WENO const& Weight_WENO,
    int Num_Stencil)
{
    f_Node = 0.0;

    for (int i = 0; i < Num_Stencil; ++i)
        f_Node += Weight_WENO[i] * f_WENO[i];
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Reconstruction2Node_WENO5_D(
    amrex::Real& f_Node,
    amrex::Real const& f_mm,
    amrex::Real const& f_m,
    amrex::Real const& f_c,
    amrex::Real const& f_p,
    amrex::Real const& f_pp,
    int ID_Node)
{
    const int Num_Stencil = 3;

    amrex::Array<amrex::Real,Num_Stencil> f_WENO;
    amrex::Array<amrex::Real,Num_Stencil> Weight_WENO;

    WENO5_Quad(f_WENO, Weight_WENO,
          f_mm, f_m, f_c, f_p, f_pp,
          ID_Node);

    Reconstruction2Node_WENO(f_Node,
                             f_WENO,
                             Weight_WENO,
                             Num_Stencil);
}



template<int NVAR>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void MatVecMultiply(amrex::Array<amrex::Real, NVAR>& Out, 
                    amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR> const& Mat, 
                    amrex::Array<amrex::Real, NVAR> const& Vec)
{
    for (int i = 0; i < NVAR; ++i) {
        Out[i] = 0.0;
        for (int j = 0; j < NVAR; ++j) {
            Out[i] += Mat[i][j] * Vec[j];
        }
    }
}

template<int NVAR>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void MatMatMultiply(amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR>& Out, 
                    amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR> const& A, 
                    amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR> const& B)
{
    for (int i = 0; i < NVAR; ++i) {
        for (int j = 0; j < NVAR; ++j) {
            Out[i][j] = 0.0;
            for (int k = 0; k < NVAR; ++k) {
                Out[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// -------------------------------------------------------------------------
// Analytic Jacobians for the 5-Equation Multiphase Model
// M = dU/dP, Minv = dP/dU
// -------------------------------------------------------------------------
template<int NVAR>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ComputeStateJacobians(amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR>& M,
                           amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR>& Minv,
                           amrex::Array<amrex::Real, NVAR> const& P,
                           Parm const& EOS_Parameter)
{
    // Zero initialize
    for (int i=0; i<NVAR; ++i) {
        for (int j=0; j<NVAR; ++j) {
            M[i][j] = 0.0;
            Minv[i][j] = 0.0;
        }
    }

    // Indices
    const int iE = 2 + AMREX_SPACEDIM; // Index for Energy/Pressure
    const int iA = 3 + AMREX_SPACEDIM; // Index for Alpha1

    amrex::Real m1 = P[0];
    amrex::Real m2 = P[1];
    amrex::Real rho = m1 + m2;
    amrex::Real p = P[iE];
    amrex::Real alpha = P[iA];

    // Extract EOS parameters (UPDATE THESE TO MATCH YOUR STRUCT)
    amrex::Real g1 = EOS_Parameter.gamma_1;
    amrex::Real g2 = EOS_Parameter.gamma_2;
    amrex::Real pinf1 = EOS_Parameter.pinf_1;
    amrex::Real pinf2 = EOS_Parameter.pinf_2;

    // Kinetic Energy K = 0.5 * (u^2 + v^2 + w^2)
    amrex::Real K = 0.0;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        K += 0.5 * P[2+d] * P[2+d];
    }

    // EOS derivatives for internal energy: rho*e(p, alpha1)
    amrex::Real xi_p = (alpha) / (g1 - 1.0) + (1.0 - alpha) / (g2 - 1.0);
    amrex::Real xi_alpha = (p + pinf1) / (g1 - 1.0) - (p + pinf2) / (g2 - 1.0);

    // ==========================================
    // Fill M (dU / dP)
    // ==========================================
    M[0][0] = 1.0; 
    M[1][1] = 1.0; 
    M[iA][iA] = 1.0;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        int iU = 2+d;
        M[iU][0] = P[iU];       // d(rho u)/dm1 = u
        M[iU][1] = P[iU];       // d(rho u)/dm2 = u
        M[iU][iU] = rho;        // d(rho u)/du = rho
        
        M[iE][iU] = rho * P[iU]; // dE/du = rho * u
    }

    M[iE][0] = K;               // dE/dm1 = K
    M[iE][1] = K;               // dE/dm2 = K
    M[iE][iE] = xi_p;           // dE/dp = d(rho*e)/dp
    M[iE][iA] = xi_alpha;       // dE/dalpha = d(rho*e)/dalpha

    // ==========================================
    // Fill Minv (dP / dU)
    // ==========================================
    Minv[0][0] = 1.0;
    Minv[1][1] = 1.0;
    Minv[iA][iA] = 1.0;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        int iU = 2+d;
        Minv[iU][0]  = -P[iU] / rho;  // du/dm1 = -u/rho
        Minv[iU][1]  = -P[iU] / rho;  // du/dm2 = -u/rho
        Minv[iU][iU] =  1.0 / rho;    // du/d(rho*u) = 1/rho
        
        Minv[iE][iU] = -P[iU] / xi_p; // dp/d(rho*u)
    }

    Minv[iE][0]  =  K / xi_p;         // dp/dm1
    Minv[iE][1]  =  K / xi_p;         // dp/dm2
    Minv[iE][iE] =  1.0 / xi_p;       // dp/dE
    Minv[iE][iA] = -xi_alpha / xi_p;  // dp/dalpha
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ConservativeCharacteristicDecomposition(
    amrex::Array<amrex::Array<amrex::Real, 4+AMREX_SPACEDIM>, 4+AMREX_SPACEDIM>& Lambda,
    amrex::Array<amrex::Array<amrex::Real, 4+AMREX_SPACEDIM>, 4+AMREX_SPACEDIM>& L_cons,
    amrex::Array<amrex::Array<amrex::Real, 4+AMREX_SPACEDIM>, 4+AMREX_SPACEDIM>& R_cons,
    amrex::Array<amrex::Real, 4+AMREX_SPACEDIM> const& U_avg,
    Parm const& EOS_Parameter,
    int ID_Dimension,
    int ID_K)
{
    constexpr int NVAR = 4+AMREX_SPACEDIM;
    
    amrex::Array<amrex::Real, NVAR> P_avg;
    Conservative2Primitive(P_avg, U_avg, EOS_Parameter);

    amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR> L_prim, R_prim;
    Primitive2CharacteristicDecomposition(Lambda, L_prim, R_prim, P_avg, EOS_Parameter, ID_Dimension, ID_K);

    amrex::Array<amrex::Array<amrex::Real, NVAR>, NVAR> M, Minv;
    ComputeStateJacobians<NVAR>(M, Minv, P_avg, EOS_Parameter);

    // Map to Conservative Space (L_cons = L_prim * Minv, R_cons = M * R_prim)
    MatMatMultiply<NVAR>(L_cons, L_prim, Minv);
    MatMatMultiply<NVAR>(R_cons, M, R_prim);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_WENO5_Quad_K(
    int i, int j, int k,
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
    amrex::Array<amrex::Real,AMREX_SPACEDIM>& Vel,
    amrex::Real& cmax,
    amrex::Array4<amrex::Real const> const& Conservative,
    Parm const& EOS_Parameter,
    int ID_Dimension,
    int ID_K,
    int ID_Riemann)
{
    using Quad = Quadrature::GaussLegendre2;
    constexpr int NVAR = 4+AMREX_SPACEDIM;

    for (int s=0; s<NVAR; s++) Flux[s]=0.0_rt;
    for (int d=0; d<AMREX_SPACEDIM; d++) Vel[d]=0.0_rt;
    cmax = 0.0_rt;

    int di=0,dj=0,ti=0,tj=0;
    if (ID_Dimension==0) { di=1; dj=0; ti=0; tj=1; }
    else if (ID_Dimension==1) { di=0; dj=1; ti=1; tj=0; }
    
    int ID_Transverse = (ID_Dimension == 0) ? 1 : 0;

    amrex::Array<amrex::Real,NVAR> U_mmm,U_mm,U_m,U_c,U_p,U_pp;
    amrex::Array<amrex::Real,NVAR> VL[5], VR[5];

    // ====================================================================
    // 1) Normal sweep: Reconstruct Conservative directly
    // ====================================================================
    for (int it=0; it<5; it++)
    {
        int offset = it-2;

        for (int s=0; s<NVAR; s++) {
            U_mmm[s]=Conservative(i-3*di+offset*ti, j-3*dj+offset*tj, k, s);
            U_mm[s] =Conservative(i-2*di+offset*ti, j-2*dj+offset*tj, k, s);
            U_m[s]  =Conservative(i-1*di+offset*ti, j-1*dj+offset*tj, k, s);
            U_c[s]  =Conservative(i+offset*ti,      j+offset*tj,      k, s);
            U_p[s]  =Conservative(i+1*di+offset*ti, j+1*dj+offset*tj, k, s);
            U_pp[s] =Conservative(i+2*di+offset*ti, j+2*dj+offset*tj, k, s);
        }

        amrex::Array<amrex::Real, NVAR> U_avg;
        for (int s=0; s<NVAR; ++s) U_avg[s] = 0.5 * (U_m[s] + U_c[s]);

        amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda, L_cons, R_cons;
        ConservativeCharacteristicDecomposition(Lambda, L_cons, R_cons, U_avg, EOS_Parameter, ID_Dimension, ID_K);

        amrex::Array<amrex::Real, NVAR> W_mmm, W_mm, W_m, W_c, W_p, W_pp;
        MatVecMultiply<NVAR>(W_mmm, L_cons, U_mmm);
        MatVecMultiply<NVAR>(W_mm,  L_cons, U_mm);
        MatVecMultiply<NVAR>(W_m,   L_cons, U_m);
        MatVecMultiply<NVAR>(W_c,   L_cons, U_c);
        MatVecMultiply<NVAR>(W_p,   L_cons, U_p);
        MatVecMultiply<NVAR>(W_pp,  L_cons, U_pp);

        amrex::Array<amrex::Real, NVAR> WL, WR;
        for (int s=0; s<NVAR; s++) {
            amrex::Real W_node;
            Reconstruction2Node_WENO5_D(W_node, W_mmm[s], W_mm[s], W_m[s], W_c[s], W_p[s], 1);  
            WL[s] = W_node;
            
            Reconstruction2Node_WENO5_D(W_node, W_mm[s], W_m[s], W_c[s], W_p[s], W_pp[s], -1); 
            WR[s] = W_node;
        }

        MatVecMultiply<NVAR>(VL[it], R_cons, WL);
        MatVecMultiply<NVAR>(VR[it], R_cons, WR);
    }

    // ====================================================================
    // 2) Transverse dimension-by-dimension to GL2 nodes
    // ====================================================================
    amrex::Array<amrex::Array<amrex::Real,NVAR>, Quad::npts> UL_q, UR_q;
    
    amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda_trans_L, L_trans_L, R_trans_L;
    amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda_trans_R, L_trans_R, R_trans_R;
    
    ConservativeCharacteristicDecomposition(Lambda_trans_L, L_trans_L, R_trans_L, VL[2], EOS_Parameter, ID_Transverse, ID_K);
    ConservativeCharacteristicDecomposition(Lambda_trans_R, L_trans_R, R_trans_R, VR[2], EOS_Parameter, ID_Transverse, ID_K);

    amrex::Array<amrex::Real, NVAR> WL_samples[5], WR_samples[5];
    for (int it=0; it<5; ++it) {
        MatVecMultiply<NVAR>(WL_samples[it], L_trans_L, VL[it]);
        MatVecMultiply<NVAR>(WR_samples[it], L_trans_R, VR[it]);
    }

    for (int q=0; q<Quad::npts; ++q)
    {
        amrex::Array<amrex::Real, NVAR> WL_q, WR_q;
        for (int s=0; s<NVAR; ++s) {
            amrex::Real W_node;
            Reconstruction2Node_WENO5_D(W_node, WL_samples[0][s], WL_samples[1][s], WL_samples[2][s], WL_samples[3][s], WL_samples[4][s], q+2); 
            WL_q[s] = W_node;
            
            Reconstruction2Node_WENO5_D(W_node, WR_samples[0][s], WR_samples[1][s], WR_samples[2][s], WR_samples[3][s], WR_samples[4][s], q+2); 
            WR_q[s] = W_node;
        }

        MatVecMultiply<NVAR>(UL_q[q], R_trans_L, WL_q);
        MatVecMultiply<NVAR>(UR_q[q], R_trans_R, WR_q);
    }

    // ====================================================================
    // 3) Riemann solve + quad integration
    // ====================================================================
    for (int q=0; q<Quad::npts; ++q)
    {
        amrex::Real w = Quad::weight(q);
        amrex::Array<amrex::Real,NVAR> Flux_q;
        amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel_q;
        amrex::Real cmax_q;

        RiemannSolver(Flux_q, Vel_q, cmax_q, UL_q[q], UR_q[q],
                      EOS_Parameter, ID_Dimension, ID_K, ID_Riemann);

        for (int s=0; s<NVAR; ++s) Flux[s] += 0.5 * w * Flux_q[s];
        for (int d=0; d<AMREX_SPACEDIM; ++d) Vel[d] += 0.5 * w * Vel_q[d];
        cmax = amrex::max(cmax, cmax_q);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_WENO5_Quad_Prim_K(
    int i, int j, int k,
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
    amrex::Array<amrex::Real,AMREX_SPACEDIM>& Vel,
    amrex::Real& cmax,
    amrex::Array4<amrex::Real const> const& Conservative,
    Parm const& EOS_Parameter,
    int ID_Dimension,
    int ID_K,
    int ID_Riemann)
{
    using Quad = Quadrature::GaussLegendre2;
    constexpr int NVAR = 4+AMREX_SPACEDIM;

    for (int s=0; s<NVAR; s++) Flux[s]=0.0_rt;
    for (int d=0; d<AMREX_SPACEDIM; d++) Vel[d]=0.0_rt;
    cmax = 0.0_rt;

    int di=0,dj=0,ti=0,tj=0;
    if (ID_Dimension==0) { di=1; dj=0; ti=0; tj=1; }
    else if (ID_Dimension==1) { di=0; dj=1; ti=1; tj=0; }
    
    int ID_Transverse = (ID_Dimension == 0) ? 1 : 0;

    amrex::Array<amrex::Real,NVAR> U_mmm,U_mm,U_m,U_c,U_p,U_pp;
    amrex::Array<amrex::Real,NVAR> VL[5], VR[5];

    // ====================================================================
    // 1) Normal sweep: Reconstruct Conservative directly
    // ====================================================================
    for (int it=0; it<5; it++)
    {
        int offset = it-2;

        for (int s=0; s<NVAR; s++) {
            U_mmm[s]=Conservative(i-3*di+offset*ti, j-3*dj+offset*tj, k, s);
            U_mm[s] =Conservative(i-2*di+offset*ti, j-2*dj+offset*tj, k, s);
            U_m[s]  =Conservative(i-1*di+offset*ti, j-1*dj+offset*tj, k, s);
            U_c[s]  =Conservative(i+offset*ti,      j+offset*tj,      k, s);
            U_p[s]  =Conservative(i+1*di+offset*ti, j+1*dj+offset*tj, k, s);
            U_pp[s] =Conservative(i+2*di+offset*ti, j+2*dj+offset*tj, k, s);
        }

        // amrex::Array<amrex::Real, NVAR> U_avg;
        // for (int s=0; s<NVAR; ++s) U_avg[s] = 0.5 * (U_m[s] + U_c[s]);

        // amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda, L_cons, R_cons;
        // ConservativeCharacteristicDecomposition(Lambda, L_cons, R_cons, U_avg, EOS_Parameter, ID_Dimension, ID_K);

        amrex::Array<amrex::Real, NVAR> W_mmm, W_mm, W_m, W_c, W_p, W_pp;
        // MatVecMultiply<NVAR>(W_mmm, L_cons, U_mmm);
        // MatVecMultiply<NVAR>(W_mm,  L_cons, U_mm);
        // MatVecMultiply<NVAR>(W_m,   L_cons, U_m);
        // MatVecMultiply<NVAR>(W_c,   L_cons, U_c);
        // MatVecMultiply<NVAR>(W_p,   L_cons, U_p);
        // MatVecMultiply<NVAR>(W_pp,  L_cons, U_pp);
        Conservative2Primitive(W_mmm, U_mmm, EOS_Parameter);
        Conservative2Primitive(W_mm, U_mm, EOS_Parameter);
        Conservative2Primitive(W_m, U_m, EOS_Parameter);
        Conservative2Primitive(W_c, U_c, EOS_Parameter);
        Conservative2Primitive(W_p, U_p, EOS_Parameter);
        Conservative2Primitive(W_pp, U_pp, EOS_Parameter);

        amrex::Array<amrex::Real, NVAR> WL, WR;
        for (int s=0; s<NVAR; s++) {
            amrex::Real W_node;
            Reconstruction2Node_WENO5_D(W_node, W_mmm[s], W_mm[s], W_m[s], W_c[s], W_p[s], 1);  
            WL[s] = W_node;
            
            Reconstruction2Node_WENO5_D(W_node, W_mm[s], W_m[s], W_c[s], W_p[s], W_pp[s], -1); 
            WR[s] = W_node;
        }

        // MatVecMultiply<NVAR>(VL[it], R_cons, WL);
        // MatVecMultiply<NVAR>(VR[it], R_cons, WR);

        VL[it] = WL;
        VR[it] = WR;

        
    }

    // ====================================================================
    // 2) Transverse dimension-by-dimension to GL2 nodes
    // ====================================================================
    amrex::Array<amrex::Array<amrex::Real,NVAR>, Quad::npts> UL_q, UR_q;
    
    // amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda_trans_L, L_trans_L, R_trans_L;
    // amrex::Array<amrex::Array<amrex::Real,NVAR>,NVAR> Lambda_trans_R, L_trans_R, R_trans_R;
    
    // ConservativeCharacteristicDecomposition(Lambda_trans_L, L_trans_L, R_trans_L, VL[2], EOS_Parameter, ID_Transverse, ID_K);
    // ConservativeCharacteristicDecomposition(Lambda_trans_R, L_trans_R, R_trans_R, VR[2], EOS_Parameter, ID_Transverse, ID_K);

    amrex::Array<amrex::Array<amrex::Real, NVAR>, 5> WL_samples, WR_samples;
    for (int it=0; it<5; ++it) {
        // MatVecMultiply<NVAR>(WL_samples[it], L_trans_L, VL[it]);
        // MatVecMultiply<NVAR>(WR_samples[it], L_trans_R, VR[it]);
        WL_samples[it] = VL[it];
        WR_samples[it] = VR[it];
    }

    for (int q=0; q<Quad::npts; ++q)
    {
        amrex::Array<amrex::Real, NSTATE> WL_q, WR_q;
        for (int s=0; s<NVAR; ++s) {
            amrex::Real W_node;
            Reconstruction2Node_WENO5_D(W_node, WL_samples[0][s], WL_samples[1][s], WL_samples[2][s], WL_samples[3][s], WL_samples[4][s], q+2); 
            WL_q[s] = W_node;
            
            Reconstruction2Node_WENO5_D(W_node, WR_samples[0][s], WR_samples[1][s], WR_samples[2][s], WR_samples[3][s], WR_samples[4][s], q+2); 
            WR_q[s] = W_node;
        }

        // MatVecMultiply<NVAR>(UL_q[q], R_trans_L, WL_q);
        // MatVecMultiply<NVAR>(UR_q[q], R_trans_R, WR_q);
        
        Primitive2Conservative(UL_q[q], WL_q, EOS_Parameter);
        Primitive2Conservative(UR_q[q], WR_q, EOS_Parameter);
    }

    // ====================================================================
    // 3) Riemann solve + quad integration
    // ====================================================================
    for (int q=0; q<Quad::npts; ++q)
    {
        amrex::Real w = Quad::weight(q);
        amrex::Array<amrex::Real,NVAR> Flux_q;
        amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel_q;
        amrex::Real cmax_q;

        RiemannSolver(Flux_q, Vel_q, cmax_q, UL_q[q], UR_q[q],
                      EOS_Parameter, ID_Dimension, ID_K, ID_Riemann);

        for (int s=0; s<NVAR; ++s) Flux[s] += 0.5 * w * Flux_q[s];
        for (int d=0; d<AMREX_SPACEDIM; ++d) Vel[d] += 0.5 * w * Vel_q[d];
        cmax = amrex::max(cmax, cmax_q);
    }
}

// void
// Compressible_PhaseField::compute_dUdt_FV (MultiFab& mf_new,
//                          MultiFab& mf_old,
//                          int lev,
//                          Real dt_lev,
//                          int ncycle,
//                          int stage,
//                          FluxRegister* fr_as_crse,
//                          FluxRegister* fr_as_fine,
//                          FluxRegister* fr_as_crse_nc,
//                          FluxRegister* fr_as_fine_nc)
// {
// BL_PROFILE("compute_dUdt_FV()");

//     MultiFab& c_max_lev = c_max[lev];
//     if (stage == 0){
//         c_max_lev.setVal(0.0);
//     }

//     auto const prob_lo = Geom(lev).ProbLoArray();
//     const Real dx = geom[lev].CellSize(0);
// #if (AMREX_SPACEDIM > 1)
//     const Real dy = geom[lev].CellSize(1);
// #else
//     const Real dy = 1.0;
// #endif
// #if (AMREX_SPACEDIM > 2)
//     const Real dz = geom[lev].CellSize(2);
// #else
//     const Real dz = 1.0;
// #endif
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};

//     Array<MultiFab,AMREX_SPACEDIM> fluxes;
// #if (NONCONSERVATIVE == true)
//     Array<MultiFab,AMREX_SPACEDIM> fluxes_nc;
// #endif
//     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
//     {
//         BoxArray ba = grids[lev];
//         ba.surroundingNodes(idim);
//         fluxes[idim].define(ba, dmap[lev], NSTATE, 0);
// #if (NONCONSERVATIVE == true)
//         fluxes_nc[idim].define(ba, dmap[lev], NC_TERMS, 0);
// #endif
//     }

//     Parm const* lparm = d_parm;
//     const int lID_COORDSYS = ID_COORDSYS;
//     amrex::Real time = t_old[lev];

//     MultiFab mf_state_ib;
//     ID_IB = 1;
//     if (ID_IB) {
//         IB::Config ib_cfg;
//         ib_cfg.enabled = true;
//         ib_cfg.phi_eps = Real(1.0e-12);
//         MultiFab& levelset = ib_levelset[lev];
//         levelset.FillBoundary(geom[lev].periodicity());
//         mf_state_ib.define(mf_old.boxArray(), mf_old.DistributionMap(), mf_old.nComp(), mf_old.nGrowVect());
//         MultiFab::Copy(mf_state_ib, mf_old, 0, 0, mf_old.nComp(), mf_old.nGrowVect());
//         FiveEqPhysicsModel fiveeq_model;
//         fiveeq_model.wall_type = static_cast<int>(FiveEqPhysicsModel::WallType::Slip);
//         IB::ApplyGhostCells(mf_state_ib, mf_old, levelset, geom[lev], ib_cfg, fiveeq_model);
//     }

// //==================================================
// // Conservative surface integral
// //==================================================
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         for (MFIter mfi(mf_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();

//             Array4<Real const> statein = ID_IB ? mf_state_ib.const_array(mfi) : mf_old.const_array(mfi);
//             // Array4<Real const> statein = mf_old.const_array(mfi);


//             Array4<Real      > dUdt      = mf_new.array(mfi);
//             Array4<Real      > c_max_new = c_max_lev.array(mfi);
// {BL_PROFILE("compute_dUdt_FV::{ computing the fluxes }");
//             for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
//                 Array4<Real> flux = fluxes[idim].array(mfi);
// #if (NONCONSERVATIVE == true)
//                 Array4<Real> fluxNC = fluxes_nc[idim].array(mfi);
// #endif
// {BL_PROFILE("compute_dUdt_FV::{ idim-direction ParallelFor }");
//                 Box bx_face = surroundingNodes(bx, idim);
//                 amrex::ParallelFor(bx_face,
//                 [=] AMREX_GPU_DEVICE (int i, int j, int k)
//                 {
// {BL_PROFILE("compute_dUdt_FV::{ computing the idim-direction fluxes }");
//                     Array<Real,NSTATE> fhat;
//                     Array<Real,NC_TERMS> fhatNC;
//                     Real cmax;
// {BL_PROFILE("compute_dUdt_FV::{ computing the idim-direction hyperbolic fluxes }");
//                     FVM_Conservative2FluxHyperbolic_K(i, j, k,
//                                                       fhat, fhatNC, cmax,
//                                                       statein, *lparm,
//                                                       dX,
//                                                       idim,
//                                                       lparm->FiniteVolume_Parm, lparm->Physics_Parm);
// }
//                     c_max_new(i,j,k,0) = std::max(cmax, c_max_new(i,j,k,0));
// #if (DIFFUSION == true)
//                     Array<Real,NSTATE> fhatD;
//                     amrex::Real dmax;
// {BL_PROFILE("compute_dUdt_FV::{ computing the idim-direction parabolic fluxes }");
//                     FDM_Conservative2FluxDiffusion_K(i, j, k,
//                                                      fhatD, dmax,
//                                                      statein,
//                                                      prob_lo,
//                                         AMREX_D_DECL(dx,dy,dz),
//                                                      time,
//                                                      *lparm,
//                                                      idim);
// }
//                     c_max_new(i,j,k,1) = std::max(dmax, c_max_new(i,j,k,1));
// #endif
//                     for (int iState = 0; iState < NSTATE; iState++){
//                         flux(i,j,k,iState) = fhat[iState];
// #if (DIFFUSION == true)
//                         flux(i,j,k,iState) += -fhatD[iState];
// #endif
//                     }
// #if (NONCONSERVATIVE == true)
//                     for (int iState = 0; iState < NC_TERMS; iState++){
//                         fluxNC(i,j,k,iState) = fhatNC[iState];
//                     }
// #endif
// }
//                 });
// }
//             } // end idim
// }

//             AMREX_D_TERM(Array4<Real const> fluxx_c = fluxes[0].const_array(mfi);,
//                          Array4<Real const> fluxy_c = fluxes[1].const_array(mfi);,
//                          Array4<Real const> fluxz_c = fluxes[2].const_array(mfi));
// {BL_PROFILE("COMPAS::compute_dUdt_FV::surface_integral()");
//             amrex::ParallelFor(bx,
//             [=] AMREX_GPU_DEVICE (int i, int j, int k)
//             {
//                 if (lID_COORDSYS == 0){
//                     for (int iState = 0; iState < NSTATE; iState++){
//                         FVM_SurfaceIntegral(i, j, k, iState,
//                                          dUdt, statein,
//                             AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
//                             AMREX_D_DECL(dx,dy,dz));
//                     }
//                 }
//                 else {
//                     FVM_SurfaceIntegral_Coord(i,j,k,
//                                            dUdt,statein,
//                               AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
//                                            lID_COORDSYS,
//                                            dX,
//                                            prob_lo);
//                 }
//             });
// }

// #ifdef USER_SOURCE_TERM
// {BL_PROFILE("COMPAS::compute_dUdt_FV::user_source_term()");
//             amrex::ParallelFor(bx,
//             [=] AMREX_GPU_DEVICE (int i, int j, int k)
//             {
//                 amrex::Real x = prob_lo[0] + (Real(i)+0.5)*dx;
// #if (AMREX_SPACEDIM > 1)
//                 amrex::Real y = prob_lo[1] + (Real(j)+0.5)*dy;
// #else
//                 amrex::Real y = 0.0;
// #endif
// #if (AMREX_SPACEDIM > 2)
//                 amrex::Real z = prob_lo[2] + (Real(k)+0.5)*dz;
// #else
//                 amrex::Real z = 0.0;
// #endif
//                 amrex::Real t = time;

//                 amrex::Array<amrex::Real, NSTATE> Residual;

//                 user_source_term(
//                     Residual,
//                     statein,
//                     i, j, k,
//                     x, y, z, t, dt_lev,
//                     AMREX_D_DECL(dx,dy,dz),
//                     *lparm);

//                 for (int iState = 0; iState < NSTATE; iState++){
//                     dUdt(i,j,k,iState) += Residual[iState];
//                 }
//             });
// }
// #endif

// #if (NONCONSERVATIVE == true)
//             AMREX_D_TERM(Array4<Real const> fluxxNC_c = fluxes_nc[0].const_array(mfi);,
//                          Array4<Real const> fluxyNC_c = fluxes_nc[1].const_array(mfi);,
//                          Array4<Real const> fluxzNC_c = fluxes_nc[2].const_array(mfi));
// {BL_PROFILE("COMPAS::compute_dUdt_FV::surface_integral_nc()");
//             amrex::ParallelFor(bx,
//             [=] AMREX_GPU_DEVICE (int i, int j, int k)
//             {
//                 FVM_SurfaceIntegral_NC(i, j, k,
//                                     dUdt, statein,
//                        AMREX_D_DECL(fluxx_c  , fluxy_c  , fluxz_c  ),
//                        AMREX_D_DECL(fluxxNC_c, fluxyNC_c, fluxzNC_c),
//                        AMREX_D_DECL(dx       , dy       , dz       ),
//                                     *lparm,
//                                     lparm->Physics_Parm);
//             });
// }
// #endif

//             if (do_reflux)
//             {
//                 for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
//                     amrex::Real dArea;
//                     if (idim == 0){
//                         dArea = dy*dz;
//                     }
//                     else if (idim == 1){
//                         dArea = dx*dz;
//                     }
//                     else {
//                         dArea = dx*dy;
//                     }

//                     Array4<Real> flux = fluxes[idim].array(mfi);
// #if (NONCONSERVATIVE == true)
//                     Array4<Real> fluxNC = fluxes_nc[idim].array(mfi);
// #endif
//                     Box bx_face = surroundingNodes(bx, idim);
//                     amrex::ParallelFor(bx_face,
//                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
//                     {
//                         for (int iState = 0; iState < NSTATE; iState++){
//                             flux(i,j,k,iState) *= dArea;
//                         }
// #if (NONCONSERVATIVE == true)
//                         for (int iState = 0; iState < NC_TERMS; iState++){
//                             fluxNC(i,j,k,iState) *= dArea;
//                         }
// #endif
//                     });
//                 } // end idim
//             }
//         } // end mfi
//     } // end omp

//     if (ID_IB) {
//         MultiFab const& levelset = ib_levelset[lev];
//         IB::ZeroResidualInSolid(mf_new, levelset, Real(0.0));
//     }

//     if (check_cfl){
// #if (!AMREX_USE_GPU)
//         Real cmax = c_max[lev].norminf(0,0,true);

//         if (AMREX_D_TERM(cmax*dt_lev > dx, ||
//                          cmax*dt_lev > dy, ||
//                          cmax*dt_lev > dz))
//         {
//             amrex::AllPrint() << "cmax = " << cmax
//                               << ", dt = " << dt_lev << " dx = " << dx << " " << dy << " " << dz << std::endl;
//             amrex::Abort("CFL violation. use smaller adv.cfl.");
//         }
// #endif
//     }

//     if (fr_as_crse) {
//         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//             const Real scale = -dt_lev;
//             fr_as_crse->CrseInit(fluxes[idim], idim, 0, 0, NSTATE, scale, FluxRegister::ADD);
// #if (NONCONSERVATIVE == true)
//             fr_as_crse_nc->CrseInit(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale, FluxRegister::ADD);
// #endif
//         }
//     }

//     if (fr_as_fine) {
//         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
//             const Real scale = dt_lev;
//             fr_as_fine->FineAdd(fluxes[idim], idim, 0, 0, NSTATE, scale);
// #if (NONCONSERVATIVE == true)
//             fr_as_fine_nc->FineAdd(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale);
// #endif
//         }
//     }
// }

#endif
#endif