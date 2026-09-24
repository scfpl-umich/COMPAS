// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ImmersedBoundaries_H_
#define ImmersedBoundaries_H_

#include <Macros.H>

using namespace amrex;

namespace IB
{

// ============================================================
// Geometry-only configuration
// ============================================================

struct Config
{
    bool enabled = true;
    Real phi_eps = Real(1.0e-12);
    bool validate_second_order_state = true;
};

// ============================================================
// Utilities
// ============================================================

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int clamp_int(int v, int lo, int hi)
{
    return amrex::max(lo, amrex::min(v, hi));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real clamp_real(Real v, Real lo, Real hi)
{
    return amrex::max(lo, amrex::min(v, hi));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_reasonable_number(Real v)
{
    return (v == v) && (amrex::Math::abs(v) < Real(1.0e150));
}

// ------------------------------------------------------------
// Identify IB ghost cell
// ------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_ib_ghost_cell(
    int i, int j, int k,
    Array4<Real const> const& phi,
    Dim3 lo, Dim3 hi,
    Real phi_eps)
{
    if (phi(i,j,k,0) <= phi_eps) {
        return false;
    }

    const int ip = clamp_int(i+1, lo.x, hi.x);
    const int im = clamp_int(i-1, lo.x, hi.x);

#if (AMREX_SPACEDIM > 1)
    const int jp = clamp_int(j+1, lo.y, hi.y);
    const int jm = clamp_int(j-1, lo.y, hi.y);
#endif
#if (AMREX_SPACEDIM > 2)
    const int kp = clamp_int(k+1, lo.z, hi.z);
    const int km = clamp_int(k-1, lo.z, hi.z);
#endif

    bool fluid = false;

    fluid |= (phi(ip,j,k,0) <= phi_eps);
    fluid |= (phi(im,j,k,0) <= phi_eps);

#if (AMREX_SPACEDIM > 1)
    fluid |= (phi(i,jp,k,0) <= phi_eps);
    fluid |= (phi(i,jm,k,0) <= phi_eps);
#endif
#if (AMREX_SPACEDIM > 2)
    fluid |= (phi(i,j,kp,0) <= phi_eps);
    fluid |= (phi(i,j,km,0) <= phi_eps);
#endif

    return fluid;
}

// ------------------------------------------------------------
// Compute interface normal and |grad(phi)|
// ------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_normal(
    int i, int j, int k,
    Array4<Real const> const& phi,
    Dim3 lo, Dim3 hi,
    GpuArray<Real,AMREX_SPACEDIM> const& dX,
    Real& nx, Real& ny, Real& nz,
    Real& gradmag)
{
    const int ip = clamp_int(i+1, lo.x, hi.x);
    const int im = clamp_int(i-1, lo.x, hi.x);

    const Real dphidx =
        (phi(ip,j,k,0) - phi(im,j,k,0)) /
        (Real(amrex::max(1,ip-im)) * dX[0]);

#if (AMREX_SPACEDIM > 1)
    const int jp = clamp_int(j+1, lo.y, hi.y);
    const int jm = clamp_int(j-1, lo.y, hi.y);

    const Real dphidy =
        (phi(i,jp,k,0) - phi(i,jm,k,0)) /
        (Real(amrex::max(1,jp-jm)) * dX[1]);
#else
    const Real dphidy = Real(0.0);
#endif

#if (AMREX_SPACEDIM > 2)
    const int kp = clamp_int(k+1, lo.z, hi.z);
    const int km = clamp_int(k-1, lo.z, hi.z);

    const Real dphidz =
        (phi(i,j,kp,0) - phi(i,j,km,0)) /
        (Real(amrex::max(1,kp-km)) * dX[2]);
#else
    const Real dphidz = Real(0.0);
#endif

    gradmag = std::sqrt(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz);

    // if (gradmag > Real(1.0e-14)) {
        nx = dphidx / gradmag;
        ny = dphidy / gradmag;
        nz = dphidz / gradmag;
    // } else {
    //     gradmag = Real(1.0);
    //     nx = Real(0.0);
    //     ny = Real(0.0);
    //     nz = Real(0.0);
    // }
}

// ------------------------------------------------------------
// 2nd-order image-state interpolation (bi-/tri-linear)
// ------------------------------------------------------------

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void sample_state_second_order(
    Array<Real,NSTATE>& Uimg,
    Array4<Real const> const& state,
    Real x, Real y, Real z,
    GpuArray<Real,AMREX_SPACEDIM> const& prob_lo,
    GpuArray<Real,AMREX_SPACEDIM> const& dX,
    Dim3 lo,
    Dim3 hi)
{
    const Real fx = (x - prob_lo[0]) / dX[0] - Real(0.5);
    int i0 = static_cast<int>(amrex::Math::floor(fx));
    Real tx = fx - Real(i0);

#if (AMREX_SPACEDIM > 1)
    const Real fy = (y - prob_lo[1]) / dX[1] - Real(0.5);
    int j0 = static_cast<int>(amrex::Math::floor(fy));
    Real ty = fy - Real(j0);
#else
    int j0 = 0;
    Real ty = Real(0.0);
#endif

#if (AMREX_SPACEDIM > 2)
    const Real fz = (z - prob_lo[2]) / dX[2] - Real(0.5);
    int k0 = static_cast<int>(amrex::Math::floor(fz));
    Real tz = fz - Real(k0);
#else
    int k0 = 0;
    Real tz = Real(0.0);
#endif

    i0 = clamp_int(i0, lo.x, hi.x-1);
    tx = clamp_real(tx, Real(0.0), Real(1.0));
    const int i1 = clamp_int(i0+1, lo.x, hi.x);

#if (AMREX_SPACEDIM > 1)
    j0 = clamp_int(j0, lo.y, hi.y-1);
    ty = clamp_real(ty, Real(0.0), Real(1.0));
    const int j1 = clamp_int(j0+1, lo.y, hi.y);
#else
    const int j1 = j0;
#endif

#if (AMREX_SPACEDIM > 2)
    k0 = clamp_int(k0, lo.z, hi.z-1);
    tz = clamp_real(tz, Real(0.0), Real(1.0));
    const int k1 = clamp_int(k0+1, lo.z, hi.z);
#else
    const int k1 = k0;
#endif

    for (int n = 0; n < NSTATE; ++n) {
#if (AMREX_SPACEDIM == 1)
        Uimg[n] = (Real(1.0)-tx) * state(i0,j0,k0,n) + tx * state(i1,j0,k0,n);
#elif (AMREX_SPACEDIM == 2)
        const Real c00 = state(i0,j0,k0,n);
        const Real c10 = state(i1,j0,k0,n);
        const Real c01 = state(i0,j1,k0,n);
        const Real c11 = state(i1,j1,k0,n);
        Uimg[n] = (Real(1.0)-tx)*(Real(1.0)-ty)*c00
                + tx*(Real(1.0)-ty)*c10
                + (Real(1.0)-tx)*ty*c01
                + tx*ty*c11;
#else
        const Real c000 = state(i0,j0,k0,n);
        const Real c100 = state(i1,j0,k0,n);
        const Real c010 = state(i0,j1,k0,n);
        const Real c110 = state(i1,j1,k0,n);
        const Real c001 = state(i0,j0,k1,n);
        const Real c101 = state(i1,j0,k1,n);
        const Real c011 = state(i0,j1,k1,n);
        const Real c111 = state(i1,j1,k1,n);

        const Real c00 = (Real(1.0)-tx)*c000 + tx*c100;
        const Real c10 = (Real(1.0)-tx)*c010 + tx*c110;
        const Real c01 = (Real(1.0)-tx)*c001 + tx*c101;
        const Real c11 = (Real(1.0)-tx)*c011 + tx*c111;

        const Real c0 = (Real(1.0)-ty)*c00 + ty*c10;
        const Real c1 = (Real(1.0)-ty)*c01 + ty*c11;

        Uimg[n] = (Real(1.0)-tz)*c0 + tz*c1;
#endif
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool second_order_stencil_all_fluid(
    Real x, Real y, Real z,
    Array4<Real const> const& phi,
    GpuArray<Real,AMREX_SPACEDIM> const& prob_lo,
    GpuArray<Real,AMREX_SPACEDIM> const& dX,
    Dim3 lo,
    Dim3 hi,
    Real phi_eps)
{
    const Real fx = (x - prob_lo[0]) / dX[0] - Real(0.5);
    int i0 = static_cast<int>(amrex::Math::floor(fx));
    i0 = clamp_int(i0, lo.x, hi.x-1);
    const int i1 = clamp_int(i0+1, lo.x, hi.x);

#if (AMREX_SPACEDIM > 1)
    const Real fy = (y - prob_lo[1]) / dX[1] - Real(0.5);
    int j0 = static_cast<int>(amrex::Math::floor(fy));
    j0 = clamp_int(j0, lo.y, hi.y-1);
    const int j1 = clamp_int(j0+1, lo.y, hi.y);
#else
    const int j0 = 0;
    const int j1 = 0;
#endif

#if (AMREX_SPACEDIM > 2)
    const Real fz = (z - prob_lo[2]) / dX[2] - Real(0.5);
    int k0 = static_cast<int>(amrex::Math::floor(fz));
    k0 = clamp_int(k0, lo.z, hi.z-1);
    const int k1 = clamp_int(k0+1, lo.z, hi.z);
#else
    const int k0 = 0;
    const int k1 = 0;
#endif

#if (AMREX_SPACEDIM == 1)
    return (phi(i0,j0,k0,0) <= phi_eps) && (phi(i1,j0,k0,0) <= phi_eps);
#elif (AMREX_SPACEDIM == 2)
    return (phi(i0,j0,k0,0) <= phi_eps) &&
           (phi(i1,j0,k0,0) <= phi_eps) &&
           (phi(i0,j1,k0,0) <= phi_eps) &&
           (phi(i1,j1,k0,0) <= phi_eps);
#else
    return (phi(i0,j0,k0,0) <= phi_eps) &&
           (phi(i1,j0,k0,0) <= phi_eps) &&
           (phi(i0,j1,k0,0) <= phi_eps) &&
           (phi(i1,j1,k0,0) <= phi_eps) &&
           (phi(i0,j0,k1,0) <= phi_eps) &&
           (phi(i1,j0,k1,0) <= phi_eps) &&
           (phi(i0,j1,k1,0) <= phi_eps) &&
           (phi(i1,j1,k1,0) <= phi_eps);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void copy_state_at_cell(
    Array<Real,NSTATE>& U,
    Array4<Real const> const& state,
    int i, int j, int k)
{
    for (int n = 0; n < NSTATE; ++n) {
        U[n] = state(i,j,k,n);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nearest_cell_from_position(
    int& ii, int& jj, int& kk,
    Real x, Real y, Real z,
    GpuArray<Real,AMREX_SPACEDIM> const& prob_lo,
    GpuArray<Real,AMREX_SPACEDIM> const& dX,
    Dim3 lo,
    Dim3 hi)
{
    const Real fx = (x - prob_lo[0]) / dX[0] - Real(0.5);
    ii = clamp_int(static_cast<int>(amrex::Math::floor(fx + Real(0.5))), lo.x, hi.x);

#if (AMREX_SPACEDIM > 1)
    const Real fy = (y - prob_lo[1]) / dX[1] - Real(0.5);
    jj = clamp_int(static_cast<int>(amrex::Math::floor(fy + Real(0.5))), lo.y, hi.y);
#else
    amrex::ignore_unused(y);
    jj = 0;
#endif

#if (AMREX_SPACEDIM > 2)
    const Real fz = (z - prob_lo[2]) / dX[2] - Real(0.5);
    kk = clamp_int(static_cast<int>(amrex::Math::floor(fz + Real(0.5))), lo.z, hi.z);
#else
    amrex::ignore_unused(z);
    kk = 0;
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sample_fluid_fallback_state(
    Array<Real,NSTATE>& Ufb,
    int i, int j, int k,
    Real xg, Real yg, Real zg,
    Real nx, Real ny, Real nz,
    Array4<Real const> const& state,
    Array4<Real const> const& phi,
    GpuArray<Real,AMREX_SPACEDIM> const& prob_lo,
    GpuArray<Real,AMREX_SPACEDIM> const& dX,
    Dim3 s_lo, Dim3 s_hi,
    Dim3 p_lo, Dim3 p_hi,
    Real phi_eps)
{
#if (AMREX_SPACEDIM == 1)
    const Real ds = dX[0];
    const int max_steps = (s_hi.x - s_lo.x + 1);
#elif (AMREX_SPACEDIM == 2)
    const Real ds = amrex::min(dX[0], dX[1]);
    const int max_steps = amrex::max(s_hi.x - s_lo.x + 1, s_hi.y - s_lo.y + 1);
#else
    const Real ds = amrex::min(dX[0], amrex::min(dX[1], dX[2]));
    const int max_steps = amrex::max(s_hi.x - s_lo.x + 1,
                          amrex::max(s_hi.y - s_lo.y + 1, s_hi.z - s_lo.z + 1));
#endif

    for (int m = 1; m <= max_steps; ++m) {
        const Real dist = (Real(m) + Real(0.5)) * ds;
        const Real xs = xg - dist * nx;
        const Real ys = yg - dist * ny;
        const Real zs = zg - dist * nz;

        int ii, jj, kk;
        nearest_cell_from_position(ii, jj, kk, xs, ys, zs, prob_lo, dX, s_lo, s_hi);

        const int ip = clamp_int(ii, p_lo.x, p_hi.x);
#if (AMREX_SPACEDIM > 1)
        const int jp = clamp_int(jj, p_lo.y, p_hi.y);
#else
        const int jp = 0;
#endif
#if (AMREX_SPACEDIM > 2)
        const int kp = clamp_int(kk, p_lo.z, p_hi.z);
#else
        const int kp = 0;
#endif

        if (phi(ip,jp,kp,0) <= phi_eps) {
            copy_state_at_cell(Ufb, state, ii, jj, kk);
            return true;
        }
    }

    // Local backup: 3x3x3 fluid neighbor with smallest phi.
    bool found = false;
    Real best_phi = Real(1.0e300);
    int bi = i, bj = j, bk = k;

    for (int ii = i-1; ii <= i+1; ++ii) {
        const int ic = clamp_int(ii, p_lo.x, p_hi.x);
#if (AMREX_SPACEDIM > 1)
        for (int jj = j-1; jj <= j+1; ++jj) {
            const int jc = clamp_int(jj, p_lo.y, p_hi.y);
#else
        for (int jj = j; jj <= j; ++jj) {
            const int jc = 0;
#endif
#if (AMREX_SPACEDIM > 2)
            for (int kk = k-1; kk <= k+1; ++kk) {
                const int kc = clamp_int(kk, p_lo.z, p_hi.z);
#else
            for (int kk = k; kk <= k; ++kk) {
                const int kc = 0;
#endif
                const Real ph = phi(ic,jc,kc,0);
                if (ph <= phi_eps && ph < best_phi) {
                    best_phi = ph;
                    bi = clamp_int(ii, s_lo.x, s_hi.x);
#if (AMREX_SPACEDIM > 1)
                    bj = clamp_int(jj, s_lo.y, s_hi.y);
#else
                    bj = 0;
#endif
#if (AMREX_SPACEDIM > 2)
                    bk = clamp_int(kk, s_lo.z, s_hi.z);
#else
                    bk = 0;
#endif
                    found = true;
                }
            }
        }
    }

    if (found) {
        copy_state_at_cell(Ufb, state, bi, bj, bk);
        return true;
    }

    return false;
}

// ============================================================
// Core IB application (policy-based)
// ============================================================

template <typename PhysicsModel>
void ApplyGhostCells(
    MultiFab& state_out,
    MultiFab const& state_src,
    MultiFab const& levelset,
    Geometry const& geom,
    Config const& cfg,
    PhysicsModel const& physics)
{
    BL_PROFILE("IB::ApplyGhostCells()");

    if (!cfg.enabled) {
        return;
    }

    AMREX_ALWAYS_ASSERT(state_out.boxArray() == state_src.boxArray());
    AMREX_ALWAYS_ASSERT(state_out.DistributionMap() == state_src.DistributionMap());
    AMREX_ALWAYS_ASSERT(levelset.boxArray() == state_src.boxArray());
    AMREX_ALWAYS_ASSERT(levelset.DistributionMap() == state_src.DistributionMap());
    AMREX_ALWAYS_ASSERT(levelset.nComp() >= 1);

    auto const prob_lo = geom.ProbLoArray();
    auto const dX = geom.CellSizeArray();
    PhysicsModel physics_model = physics;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state_out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto sout = state_out.array(mfi);
        auto ssrc = state_src.const_array(mfi);
        auto phi  = levelset.const_array(mfi);

        const Dim3 s_lo = lbound(state_src[mfi].box());
        const Dim3 s_hi = ubound(state_src[mfi].box());
        const Dim3 p_lo = lbound(levelset[mfi].box());
        const Dim3 p_hi = ubound(levelset[mfi].box());

        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            if (!is_ib_ghost_cell(i,j,k,phi,p_lo,p_hi,cfg.phi_eps)) {
                return;
            }

            Real nx, ny, nz, gradmag;
            compute_normal(i,j,k,phi,p_lo,p_hi,dX,nx,ny,nz,gradmag);

            const Real xg = prob_lo[0] + (Real(i) + Real(0.5)) * dX[0];
#if (AMREX_SPACEDIM > 1)
            const Real yg = prob_lo[1] + (Real(j) + Real(0.5)) * dX[1];
#else
            const Real yg = Real(0.0);
#endif
#if (AMREX_SPACEDIM > 2)
            const Real zg = prob_lo[2] + (Real(k) + Real(0.5)) * dX[2];
#else
            const Real zg = Real(0.0);
#endif

            // Geometry-consistent projection distance.
            const Real d = amrex::max(phi(i,j,k,0) / amrex::max(gradmag, Real(1.0e-14)), Real(0.0));
            const Real x_img = xg - Real(2.0) * d * nx;
            const Real y_img = yg - Real(2.0) * d * ny;
            const Real z_img = zg - Real(2.0) * d * nz;

            Array<Real,NSTATE> Uimg;
            bool have_img = false;

            const bool fluid_stencil = second_order_stencil_all_fluid(
                x_img, y_img, z_img, phi, prob_lo, dX, p_lo, p_hi, cfg.phi_eps);

            if (fluid_stencil) {
                sample_state_second_order(Uimg, ssrc, x_img, y_img, z_img, prob_lo, dX, s_lo, s_hi);
                // have_img = (!cfg.validate_second_order_state) || physics_model.is_admissible(Uimg);
                have_img = true;
            }

            if (!have_img) {
                Array<Real,NSTATE> Ufb;
                const bool fb_ok = sample_fluid_fallback_state(
                    Ufb,
                    i, j, k,
                    xg, yg, zg,
                    nx, ny, nz,
                    ssrc, phi,
                    prob_lo, dX,
                    s_lo, s_hi,
                    p_lo, p_hi,
                    cfg.phi_eps);
                if (fb_ok) {
                    Uimg = Ufb;
                    have_img = true;
                }
            }

            if (!have_img) {
                // No fluid sample found: preserve existing value.
                return;
            }

            Array<Real,NSTATE> Ughost;
            physics_model.apply_wall(Ughost, Uimg, nx, ny, nz);

            for (int n = 0; n < NSTATE; ++n) {
                sout(i,j,k,n) = Ughost[n];
            }
        });
    }

    state_out.FillBoundary(geom.periodicity());
}

inline void ZeroResidualInSolid(
    MultiFab& residual,
    MultiFab const& levelset,
    Real phi_eps)
{
    BL_PROFILE("IB::ZeroResidualInSolid()");

    AMREX_ALWAYS_ASSERT(levelset.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(levelset.boxArray() == residual.boxArray());
    AMREX_ALWAYS_ASSERT(levelset.DistributionMap() == residual.DistributionMap());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(residual, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto res = residual.array(mfi);
        auto phi = levelset.const_array(mfi);
        const int ncomp = residual.nComp();

        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i,int j,int k)
        {
            if (phi(i,j,k,0) > phi_eps)
            {
                for (int n=0; n<ncomp; ++n) {
                    res(i,j,k,n) = Real(0.0);
                }
            }
        });
    }
}

} // namespace IB

struct FiveEqPhysicsModel
{
    enum class WallType : int { Slip = 0, NoSlip = 1 };

    Real rho_floor = Real(1.0e-14);
    Real mass_floor = Real(1.0e-14);
    Real eint_floor = Real(1.0e-14);

    int wall_type = static_cast<int>(WallType::Slip);

    int mass1_comp = INDEX_Mass1;
    int mass2_comp = INDEX_Mass2;
    int momx_comp = INDEX_MomentumX;
#if (AMREX_SPACEDIM > 1)
    int momy_comp = INDEX_MomentumY;
#else
    int momy_comp = -1;
#endif
#if (AMREX_SPACEDIM > 2)
    int momz_comp = INDEX_MomentumZ;
#else
    int momz_comp = -1;
#endif
    int E_comp = INDEX_Energy;
    int alpha1_comp = INDEX_VolumeFraction1;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool is_admissible(Array<Real,NSTATE> const& U) const
    {
        for (int n = 0; n < NSTATE; ++n) {
            if (!IB::is_reasonable_number(U[n])) {
                return false;
            }
        }

        const Real m1 = U[mass1_comp];
        const Real m2 = U[mass2_comp];
        const Real rho = m1 + m2;
        // if (m1 <= mass_floor || m2 <= mass_floor || rho <= rho_floor) {
        //     return false;
        // }

        if (alpha1_comp >= 0 && alpha1_comp < NSTATE) {
            const Real a1 = U[alpha1_comp];
            // if (a1 < Real(-1.0e-6) || a1 > Real(1.0) + Real(1.0e-6)) {
            //     return false;
            // }
        }

        if (E_comp >= 0 && E_comp < NSTATE) {
            const Real invrho = Real(1.0) / rho;
            const Real ux = U[momx_comp] * invrho;
            Real uy = Real(0.0);
            Real uz = Real(0.0);
#if (AMREX_SPACEDIM > 1)
            uy = U[momy_comp] * invrho;
#endif
#if (AMREX_SPACEDIM > 2)
            uz = U[momz_comp] * invrho;
#endif
            const Real ke = Real(0.5) * rho * (ux*ux + uy*uy + uz*uz);
            // if (U[E_comp] <= ke) {
            //     return false;
            // }
        }

        return true;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void apply_wall(
        Array<Real,NSTATE>& Ughost,
        Array<Real,NSTATE> const& Uimg,
        Real nx, Real ny, Real nz) const
    {
        for (int n = 0; n < NSTATE; ++n) {
            Ughost[n] = Uimg[n];
        }

        const Real m1 = Uimg[mass1_comp];
        const Real m2 = Uimg[mass2_comp];
        const Real rho = amrex::max(m1 + m2, rho_floor);
        const Real invrho = Real(1.0) / rho;

        Real ux = Uimg[momx_comp] * invrho;
        Real uy = Real(0.0);
        Real uz = Real(0.0);
#if (AMREX_SPACEDIM > 1)
        uy = Uimg[momy_comp] * invrho;
#endif
#if (AMREX_SPACEDIM > 2)
        uz = Uimg[momz_comp] * invrho;
#endif

        Real ugx = ux;
        Real ugy = uy;
        Real ugz = uz;

        if (wall_type == static_cast<int>(WallType::Slip)) {
            const Real un = ux*nx + uy*ny + uz*nz;
            ugx = ux - Real(2.0) * un * nx;
            ugy = uy - Real(2.0) * un * ny;
            ugz = uz - Real(2.0) * un * nz;
        } else {
            ugx = -ux;
            ugy = -uy;
            ugz = -uz;
        }

        Ughost[momx_comp] = rho * ugx;
#if (AMREX_SPACEDIM > 1)
        Ughost[momy_comp] = rho * ugy;
#endif
#if (AMREX_SPACEDIM > 2)
        Ughost[momz_comp] = rho * ugz;
#endif

        if (E_comp >= 0 && E_comp < NSTATE) {
            const Real ke_img = Real(0.5) * rho * (ux*ux + uy*uy + uz*uz);
            const Real eint = amrex::max(Uimg[E_comp] - ke_img, eint_floor);
            const Real ke_ghost = Real(0.5) * rho * (ugx*ugx + ugy*ugy + ugz*ugz);
            Ughost[E_comp] = eint + ke_ghost;
        }

        if (alpha1_comp >= 0 && alpha1_comp < NSTATE) {
            Ughost[alpha1_comp] = Uimg[alpha1_comp];
        }
        Ughost[mass1_comp] = m1;
        Ughost[mass2_comp] = m2;
    }
};

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

