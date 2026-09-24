// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_User(int i, int j, int k, int n,
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
    //qbc = BC_Function or Normal_dot_Gradqbc = BC_Function;
    //BC_Function can be a function of not only xBC, yBC, zBC, and time but also q and their derivatives at the boundaries;

    if (i < ilo) // bc at xlo
    {
        BC_Function  = 0.0;
    }
    else if (i > ihi) // bc at xhi
    {
        BC_Function  = 0.0;
    }
#if (AMREX_SPACEDIM > 1)
    if (j < jlo) // bc at ylo
    {
        BC_Function  = 0.0;
    }
    else if (j > jhi) // bc at yhi
    {
        BC_Function  = 0.0;
    }
#endif
#if (AMREX_SPACEDIM > 2)
    if (k < klo) // bc at zlo
    {
        BC_Function  = 0.0;
    }
    else if (k > khi) // bc at zhi
    {
        BC_Function  = 0.0;
    }
#endif
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_Robin(int i, int j, int k, int n,
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
    //BC_Dirichlet*qbc + BC_Neumann*Normal_dot_Gradqbc = BC_Function;
    //BC_Dirichlet, BC_Neumann, and BC_Function can be functions of not only xBC, yBC, zBC, and time but also q and their derivatives at the boundaries;

    if (i < ilo) // bc at xlo
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
    else if (i > ihi) // bc at xhi
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
#if (AMREX_SPACEDIM > 1)
    if (j < jlo) // bc at ylo
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
    else if (j > jhi) // bc at yhi
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
#endif
#if (AMREX_SPACEDIM > 2)
    if (k < klo) // bc at zlo
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
    else if (k > khi) // bc at zhi
    {
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.0;
    }
#endif
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_User(const IntVect& iv,
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
    //User defined ghost cell filling

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
    Real x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
    Real y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
    Real y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
    Real z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
    Real z = 0.0;
#endif

    for (int n = dcomp; n < numcomp+dcomp; ++n)
    {
        const BCRec& bc = bcr[bcomp+n-dcomp];

        Real xBC, yBC, zBC;
        if (i < ilo)
        {
            if (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::ext_dir_cc) {
                xBC = problo[0];
                yBC = y;
                zBC = z;
                //q(i,j,k,n) = ...;
            }
        }
        else if (i > ihi)
        {
            if (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::ext_dir_cc) {
                xBC = probhi[0];
                yBC = y;
                zBC = z;
                //q(i,j,k,n) = ...;
            }
        }
#if (AMREX_SPACEDIM > 1)
        if (j < jlo)
        {
            if (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::ext_dir_cc) {
                xBC = x;
                yBC = problo[1];
                zBC = z;
                //q(i,j,k,n) = ...;
            }
        }
        else if (j > jhi)
        {
            if (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::ext_dir_cc) {
                xBC = x;
                yBC = probhi[1];
                zBC = z;
                //q(i,j,k,n) = ...;
            }
        }
#endif
#if (AMREX_SPACEDIM > 2)
        if (k < klo)
        {
            if (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::ext_dir_cc) {
                xBC = x;
                yBC = y;
                zBC = problo[2];
                //q(i,j,k,n) = ...;
            }
        }
        else if (k > khi)
        {
            if (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::ext_dir_cc) {
                xBC = x;
                yBC = y;
                zBC = probhi[2];
                //q(i,j,k,n) = ...;
            }
        }
#endif
    }
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void BoundaryFill_User(Box const& bx,
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
    //User defined ghost cell filling

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

    for (int n = 0; n < numcomp; ++n)
    {
        Array4<Real> q(qn,n);
        BCRec const& bc = bcn[n];

        Real x, y, z;
        Real xBC, yBC, zBC;
        if (lo.x < ilo) {
           const int imin = lo.x;
           const int imax = ilo-1;
           if (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::ext_dir_cc) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = problo[0];
                    yBC = y;
                    zBC = z;
                    //q(i,j,k) = ...;
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
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = probhi[0];
                    yBC = y;
                    zBC = z;
                    //q(i,j,k) = ...;
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
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = x;
                    yBC = problo[1];
                    zBC = z;
                    //q(i,j,k) = ...;
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
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = x;
                    yBC = probhi[1];
                    zBC = z;
                    //q(i,j,k) = ...;
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
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = x;
                    yBC = y;
                    zBC = problo[2];
                    //q(i,j,k) = ...;
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
                    x = problo[0] + ( Real(i) + 0.5 )*dX[0];
#if (AMREX_SPACEDIM > 1)
                    y = problo[1] + ( Real(j) + 0.5 )*dX[1];
#else
                    y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                    z = problo[2] + ( Real(k) + 0.5 )*dX[2];
#else
                    z = 0.0;
#endif
                    xBC = x;
                    yBC = y;
                    zBC = probhi[2];
                    //q(i,j,k) = ...;
                }}}
            }
        }
#endif
    }
}