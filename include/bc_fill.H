// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.


#ifndef BCFILL_H
#define BCFILL_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include "ProblemICBC.H"

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction(int i, int j, int k, int n,
                      Real& BC_Dirichlet,
                      Real& BC_Neumann,
                      Real& BC_Function,
                      Real const& xBC,
                      Real const& yBC,
                      Real const& zBC,
                      Real const& time,
                      Array4<Real const> const& q,
         AMREX_D_DECL(int ilo, int jlo, int klo),
         AMREX_D_DECL(int ihi, int jhi, int khi),
         AMREX_D_DECL(Real const& dx,
                      Real const& dy,
                      Real const& dz),
                      Parm const* parm)
{
#if defined(USER_BOUNDARY_FUNC_ROBIN)
    BoundaryFunction_Robin(i, j, k, n,
                           BC_Dirichlet, BC_Neumann, BC_Function,
                           xBC, yBC, zBC, time,
                           q,
              AMREX_D_DECL(ilo, jlo, klo),
              AMREX_D_DECL(ihi, jhi, khi),
              AMREX_D_DECL( dx,  dy,  dz),
                           parm);
#elif defined(USER_BOUNDARY_FUNC_DIRICHLET)
    BC_Dirichlet = 1.0;
    BC_Neumann   = 0.0;
    BoundaryFunction_User(i, j, k, n,
                          BC_Function,
                          xBC, yBC, zBC, time,
                          q,
             AMREX_D_DECL(ilo, jlo, klo),
             AMREX_D_DECL(ihi, jhi, khi),
             AMREX_D_DECL( dx,  dy,  dz),
                          parm);
#elif defined(USER_BOUNDARY_FUNC_NEUMANN)
    BC_Dirichlet = 0.0;
    BC_Neumann   = 1.0;
    BoundaryFunction_User(i, j, k, n,
                          BC_Function,
                          xBC, yBC, zBC, time,
                          q,
             AMREX_D_DECL(ilo, jlo, klo),
             AMREX_D_DECL(ihi, jhi, khi),
             AMREX_D_DECL( dx,  dy,  dz),
                          parm);
#else
    Abort("User boundary function has not been defined! Define it in ProblemICBC.H\n");
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryCondition_Robin(int i, int j, int k, int n,
                             Real& qGhost,
                             Real const& qInterior,
                             Real const& Delta_X,
                             Real const& xBC,
                             Real const& yBC,
                             Real const& zBC,
                             Real const& time,
                             Array4<Real const> const& q,
                AMREX_D_DECL(int ilo, int jlo, int klo),
                AMREX_D_DECL(int ihi, int jhi, int khi),
                AMREX_D_DECL(Real const& dx,
                             Real const& dy,
                             Real const& dz),
                             Parm const* parm)
{
    Real BC_Dirichlet, BC_Neumann, BC_Function;
    BoundaryFunction(i, j, k, n,
                     BC_Dirichlet, BC_Neumann, BC_Function,
                     xBC, yBC, zBC, time,
                     q,
        AMREX_D_DECL(ilo, jlo, klo),
        AMREX_D_DECL(ihi, jhi, khi),
        AMREX_D_DECL( dx,  dy,  dz),
                     parm);

    Real bc_dirichlet = BC_Dirichlet/2.0;
    Real bc_neumann   = BC_Neumann/Delta_X;
    Real bc_function  = BC_Function;
    qGhost = (bc_function - (bc_dirichlet - bc_neumann)*qInterior)/(bc_dirichlet + bc_neumann);
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_Robin(int i, int j, int k, int n,
                        Array4<Real> const& q,
           AMREX_D_DECL(int ilo, int jlo, int klo),
           AMREX_D_DECL(int ihi, int jhi, int khi),
           AMREX_D_DECL(Real const& xlo, Real const& ylo, Real const& zlo),
           AMREX_D_DECL(Real const& xhi, Real const& yhi, Real const& zhi),
           AMREX_D_DECL(Real const&  dx, Real const&  dy, Real const&  dz),
                        Real const& time,
                        const BCRec& bc,
                        const int orig_comp,
                        Parm const* parm)
{
    //Normal_dot_Gradqbc = (q(i,j,k,n) - qN)/Delta_X;
    //qbc                = (q(i,j,k,n) + qN)/2.0;
    //BC_Dirichlet*qbc + BC_Neumann*Normal_dot_Gradqbc = BC_Function;

    Real x = xlo + ( Real(i) + 0.5 )*dx;
#if (AMREX_SPACEDIM > 1)
    Real y = ylo + ( Real(j) + 0.5 )*dy;
#else
    Real y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
    Real z = zlo + ( Real(k) + 0.5 )*dz;
#else
    Real z = 0.0;
#endif

    int iInterior, jInterior, kInterior;
    Real Delta_X;
    Real xBC, yBC, zBC;
    Real qInterior, qGhost;
    if (i < ilo)
    {
        if (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::ext_dir_cc) {
            iInterior = 2*ilo - i - 1;
            jInterior = j;
            kInterior = k;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(iInterior - i)*dx;

            xBC = xlo;
            yBC = y;
            zBC = z;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
    else if (i > ihi)
    {
        if (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::ext_dir_cc) {
            iInterior = 2*ihi - i + 1;
            jInterior = j;
            kInterior = k;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(i - iInterior)*dx;

            xBC = xhi;
            yBC = y;
            zBC = z;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
#if (AMREX_SPACEDIM > 1)
    if (j < jlo)
    {
        if (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::ext_dir_cc) {
            iInterior = i;
            jInterior = 2*jlo - j - 1;
            kInterior = k;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(jInterior - j)*dy;

            xBC = x;
            yBC = ylo;
            zBC = z;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
    else if (j > jhi)
    {
        if (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::ext_dir_cc) {
            iInterior = i;
            jInterior = 2*jhi - j + 1;
            kInterior = k;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(j - jInterior)*dy;

            xBC = x;
            yBC = yhi;
            zBC = z;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
#endif
#if (AMREX_SPACEDIM > 2)
    if (k < klo)
    {
        if (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::ext_dir_cc) {
            iInterior = i;
            jInterior = j;
            kInterior = 2*klo - k - 1;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(kInterior - k)*dz;

            xBC = x;
            yBC = y;
            zBC = zlo;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
    else if (k > khi)
    {
        if (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::ext_dir_cc) {
            iInterior = i;
            jInterior = j;
            kInterior = 2*khi - k + 1;
            qInterior = q(iInterior, jInterior, kInterior, n);
            Delta_X = Real(k - kInterior)*dz;

            xBC = x;
            yBC = y;
            zBC = zhi;
            BoundaryCondition_Robin(i, j, k, n,
                                    qGhost,
                                    qInterior, Delta_X,
                                    xBC, yBC, zBC, time,
                                    q,
                       AMREX_D_DECL(ilo, jlo, klo),
                       AMREX_D_DECL(ihi, jhi, khi),
                       AMREX_D_DECL( dx,  dy,  dz),
                                    parm);
            q(i,j,k,n) = qGhost;
        }
    }
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_Robin(const IntVect& iv,
                        Array4<Real> const& q,
                        const int dcomp,
                        const int numcomp,
                        GeometryData const& geom,
                        const Real time,
                        const BCRec* bcr,
                        const int bcomp,
                        const int orig_comp,
                        Parm const* parm)
{
    //Normal_dot_Gradqbc = (q(i,j,k,n) - qN)/Delta_X;
    //qbc                = (q(i,j,k,n) + qN)/2.0;
    //BC_Dirichlet*qbc + BC_Neumann*Normal_dot_Gradqbc = BC_Function;

    const int i = iv[0];
#if (AMREX_SPACEDIM > 1)
    const int j = iv[1];
#else
    const int j = 0;
#endif
#if (AMREX_SPACEDIM > 2)
    const int k = iv[2];
#else
    const int k = 0;
#endif

    const auto& domain_box = geom.domain;
    const auto& domain_lo = domain_box.loVect();
    const auto& domain_hi = domain_box.hiVect();
    const int ilo = domain_lo[0];
    const int ihi = domain_hi[0];
#if (AMREX_SPACEDIM > 1)
    const int jlo = domain_lo[1];
    const int jhi = domain_hi[1];
#endif
#if (AMREX_SPACEDIM > 2)
    const int klo = domain_lo[2];
    const int khi = domain_hi[2];
#endif

    const auto dX = geom.CellSize();
    const auto problo = geom.ProbLo();
    const auto probhi = geom.ProbHi();

    for (int n = dcomp; n < numcomp+dcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n-dcomp];

        BoundaryFill_Robin(i, j, k, n,
                           q,
              AMREX_D_DECL(ilo, jlo, klo),
              AMREX_D_DECL(ihi, jhi, khi),
              AMREX_D_DECL(problo[0], problo[1], problo[2]),
              AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
              AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                           time,
                           bc,
                           orig_comp,
                           parm);
    }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_Robin(Box const& bx,
                        Array4<Real> const& qn,
                        const int dcomp,
                        const int numcomp,
                        GeometryData const& geom,
                        const Real time,
                        const BCRec* bcn,
                        const int bcomp,
                        const int orig_comp,
                        Parm const* parm)
{
    //Normal_dot_Gradqbc = (q(i,j,k) - qN)/Delta_x;
    //qbc                = (q(i,j,k) + qN)/2.0;
    //BC_Dirichlet*qbc + BC_Neumann*Normal_dot_Gradqbc = BC_Function;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Box qbx(qn);
    const auto qlo = amrex::lbound(qbx);
    const auto qhi = amrex::ubound(qbx);

    const auto domain = geom.domain;
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    const int ilo = domlo.x;
    const int ihi = domhi.x;
    const int is = std::max(qlo.x, ilo);
    const int ie = std::min(qhi.x, ihi);

#if (AMREX_SPACEDIM > 1)
    const int jlo = domlo.y;
    const int jhi = domhi.y;
    const int js = std::max(qlo.y, jlo);
    const int je = std::min(qhi.y, jhi);
#endif

#if (AMREX_SPACEDIM > 2)
    const int klo = domlo.z;
    const int khi = domhi.z;
    const int ks = std::max(qlo.z, klo);
    const int ke = std::min(qhi.z, khi);
#endif

    const auto dX = geom.CellSize();
    const auto problo = geom.ProbLo();
    const auto probhi = geom.ProbHi();

    for (int n = dcomp; n < numcomp+dcomp; ++n)
    {
        BCRec const& bc = bcn[n-dcomp];

        if (lo.x < ilo) {
           const int imin = lo.x;
           const int imax = ilo-1;
           if (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::ext_dir_cc) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
               }}}
           }
        }

        if (hi.x > ihi) {
            const int imin = ihi+1;
            const int imax = hi.x;
            if (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::ext_dir_cc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
                }}}
            }
        }

#if (AMREX_SPACEDIM > 1)
        if (lo.y < jlo) {
            const int jmin = lo.y;
            const int jmax = jlo-1;
            if (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::ext_dir_cc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
                }}}
            }
        }

        if (hi.y > jhi) {
            const int jmin = jhi+1;
            const int jmax = hi.y;
            if (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::ext_dir_cc) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
                }}}
            }
        }
#endif

#if (AMREX_SPACEDIM > 2)
        if (lo.z < klo) {
            const int kmin = lo.z;
            const int kmax = klo-1;
            if (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::ext_dir_cc) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
                }}}
            }
        }

        if (hi.z > khi) {
            const int kmin = khi+1;
            const int kmax = hi.z;
            if (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::ext_dir_cc) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    BoundaryFill_Robin(i, j, k, n,
                                       qn,
                          AMREX_D_DECL(ilo, jlo, klo),
                          AMREX_D_DECL(ihi, jhi, khi),
                          AMREX_D_DECL(problo[0], problo[1], problo[2]),
                          AMREX_D_DECL(probhi[0], probhi[1], probhi[2]),
                          AMREX_D_DECL(    dX[0],     dX[1],     dX[2]),
                                       time,
                                       bc,
                                       orig_comp,
                                       parm);
                }}}
            }
        }
#endif
    }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill(const IntVect& iv,
                  Array4<Real> const& q,
                  const int dcomp,
                  const int numcomp,
                  GeometryData const& geom,
                  const Real time,
                  const BCRec* bcr,
                  const int bcomp,
                  const int orig_comp,
                  Parm const* parm)
{
#ifdef USER_BOUNDARY_FUNC
    BoundaryFill_User(iv, q, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp, parm);
#else
    BoundaryFill_Robin(iv, q, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp, parm);
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill(Box const& bx,
                  Array4<Real> const& qn,
                  const int dcomp,
                  const int numcomp,
                  GeometryData const& geom,
                  const Real time,
                  const BCRec* bcn,
                  const int bcomp,
                  const int orig_comp,
                  Parm const* parm)
{
#ifdef USER_BOUNDARY_FUNC
    // The case defines the per-cell BoundaryFill_User; apply it to each ghost cell of bx
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        const IntVect iv(AMREX_D_DECL(i, j, k));
        if (!geom.domain.contains(iv)) {
            BoundaryFill_User(iv, qn, dcomp, numcomp, geom, time, bcn, bcomp, orig_comp, parm);
        }
    }}}
#else
    BoundaryFill_Robin(bx, qn, dcomp, numcomp, geom, time, bcn, bcomp, orig_comp, parm);
#endif
}


namespace {

struct AmrCoreGPUFill
{
    Parm const* parm;

    AMREX_GPU_HOST_DEVICE
    void operator() (const IntVect& iv,
                     Array4<Real> const& data,
                     const int dcomp,
                     const int numcomp,
                     GeometryData const& geom,
                     const Real time,
                     const BCRec* bcr,
                     const int bcomp,
                     const int orig_comp) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
        BoundaryFill(iv, data, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp, parm);
    }
};


struct AmrCoreCPUFill
{
    inline static Parm const* parm = nullptr;

    AMREX_GPU_HOST_DEVICE
    static void Boundary_Fill(Box const& bx,
                              Array4<Real> const& dest,
                              const int dcomp,
                              const int numcomp,
                              GeometryData const& geom,
                              const Real time,
                              const BCRec* bcr,
                              const int bcomp,
                              const int orig_comp)
    {
        // do something for external Dirichlet (BCType::ext_dir) if there are
        BoundaryFill(bx, dest, dcomp, numcomp, geom, time, bcr, bcomp, orig_comp, parm);
    }
};


}

#endif