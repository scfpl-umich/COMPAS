// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>


#include "ProblemICBC.H"
#include "Macros.H"


template <class Quad = Quadrature::GaussLegendre4>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_dof_ca(amrex::Box const& bx,
           amrex::Array4<amrex::Real> const& dof,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx_FinestLevel,
           int const& nvar,
           Parm const& parm,           
           bool use_quadrature = true,
           int ID_IB = 0)
{
    using namespace amrex;

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    constexpr int dim = AMREX_SPACEDIM;
    const Real inv_norm = 1.0_rt / static_cast<Real>(1 << dim); // 1/(2^dim)

    constexpr int npts = Quad::npts;

    for (int k = lo.z; k <= hi.z; ++k)
    {
#if (AMREX_SPACEDIM == 3)
        Real zc = prob_lo[2] + (0.5_rt + static_cast<Real>(k)) * dx[2];
#else
        Real zc = 0.0_rt;
#endif

        for (int j = lo.y; j <= hi.y; ++j)
        {
#if (AMREX_SPACEDIM >= 2)
            Real yc = prob_lo[1] + (0.5_rt + static_cast<Real>(j)) * dx[1];
#else
            Real yc = 0.0_rt;
#endif

            for (int i = lo.x; i <= hi.x; ++i)
            {
                Real xc = prob_lo[0] + (0.5_rt + static_cast<Real>(i)) * dx[0];

                // ---------------------------------------------------------
                // Center-sample fallback
                // ---------------------------------------------------------
                if (!use_quadrature)
                {
                    if (ID_IB)
                    {
#ifdef USE_IB
                        Real phi;
                        immersed_boundary_levelset(phi, xc, yc, zc,
                                                   AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                                   AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                                   parm);

                        dof(i,j,k,0) = phi;
#endif
                    }
                    else
                    {
                        amrex::Array<Real, NSTATE> Conservative;
                        bool sharp;

                        initial_condition(Conservative, sharp,
                                          xc, yc, zc,
                                          AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                          AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                          parm);

                        for (int v = 0; v < nvar; ++v)
                        {
                            dof(i,j,k,v) = Conservative[v];
                        }
                    }

                    continue;
                }

                // ---------------------------------------------------------
                // Quadrature projection
                // ---------------------------------------------------------

                if (ID_IB)
                {
#ifdef USE_IB
                    Real accum = 0.0_rt;

#if (AMREX_SPACEDIM == 3)
                    for (int qz = 0; qz < npts; ++qz)
                    {
                        Real rz = Quad::node(qz);
                        Real wz = Quad::weight(qz);
                        Real zq = zc + 0.5_rt * dx[2] * rz;

                        for (int qy = 0; qy < npts; ++qy)
                        {
                            Real ry = Quad::node(qy);
                            Real wy = Quad::weight(qy);
                            Real yq = yc + 0.5_rt * dx[1] * ry;

                            for (int qx = 0; qx < npts; ++qx)
                            {
                                Real rx = Quad::node(qx);
                                Real wx = Quad::weight(qx);
                                Real xq = xc + 0.5_rt * dx[0] * rx;

                                Real phi;
                                immersed_boundary_levelset(phi, xq, yq, zq,
                                                           AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                                           AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                                           parm);

                                accum += wx * wy * wz * phi;
                            }
                        }
                    }
#elif (AMREX_SPACEDIM == 2)
                    for (int qy = 0; qy < npts; ++qy)
                    {
                        Real ry = Quad::node(qy);
                        Real wy = Quad::weight(qy);
                        Real yq = yc + 0.5_rt * dx[1] * ry;

                        for (int qx = 0; qx < npts; ++qx)
                        {
                            Real rx = Quad::node(qx);
                            Real wx = Quad::weight(qx);
                            Real xq = xc + 0.5_rt * dx[0] * rx;

                            Real phi;
                            immersed_boundary_levelset(phi, xq, yq, 0.0_rt,
                                                       AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                                       AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                                       parm);

                            accum += wx * wy * phi;
                        }
                    }
#else
                    for (int qx = 0; qx < npts; ++qx)
                    {
                        Real rx = Quad::node(qx);
                        Real wx = Quad::weight(qx);
                        Real xq = xc + 0.5_rt * dx[0] * rx;

                        Real phi;
                        immersed_boundary_levelset(phi, xq, 0.0_rt, 0.0_rt,
                                                   AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                                   AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                                   parm);

                        accum += wx * phi;
                    }
#endif

                    dof(i,j,k,0) = accum * inv_norm;
#endif
                }
                else
                {
                    amrex::Array<Real, NSTATE> accum;
                    for (int v = 0; v < nvar; ++v)
                        accum[v] = 0.0_rt;

                    amrex::Array<Real, NSTATE> Conservative;
                    bool sharp;

#if (AMREX_SPACEDIM == 3)
                    for (int qz = 0; qz < npts; ++qz)
                    {
                        Real rz = Quad::node(qz);
                        Real wz = Quad::weight(qz);
                        Real zq = zc + 0.5_rt * dx[2] * rz;

                        for (int qy = 0; qy < npts; ++qy)
                        {
                            Real ry = Quad::node(qy);
                            Real wy = Quad::weight(qy);
                            Real yq = yc + 0.5_rt * dx[1] * ry;

                            for (int qx = 0; qx < npts; ++qx)
                            {
                                Real rx = Quad::node(qx);
                                Real wx = Quad::weight(qx);
                                Real xq = xc + 0.5_rt * dx[0] * rx;

                                initial_condition(Conservative, sharp,
                                                  xq, yq, zq,
                                                  AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                                  AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                                  parm);

                                for (int v = 0; v < nvar; ++v)
                                    accum[v] += wx * wy * wz * Conservative[v];
                            }
                        }
                    }
#elif (AMREX_SPACEDIM == 2)
                    for (int qy = 0; qy < npts; ++qy)
                    {
                        Real ry = Quad::node(qy);
                        Real wy = Quad::weight(qy);
                        Real yq = yc + 0.5_rt * dx[1] * ry;

                        for (int qx = 0; qx < npts; ++qx)
                        {
                            Real rx = Quad::node(qx);
                            Real wx = Quad::weight(qx);
                            Real xq = xc + 0.5_rt * dx[0] * rx;

                            initial_condition(Conservative, sharp,
                                              xq, yq, 0.0_rt,
                                              AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                              AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                              parm);

                            for (int v = 0; v < nvar; ++v)
                                accum[v] += wx * wy * Conservative[v];
                        }
                    }
#else
                    for (int qx = 0; qx < npts; ++qx)
                    {
                        Real rx = Quad::node(qx);
                        Real wx = Quad::weight(qx);
                        Real xq = xc + 0.5_rt * dx[0] * rx;

                        initial_condition(Conservative, sharp,
                                          xq, 0.0_rt, 0.0_rt,
                                          AMREX_D_DECL(dx[0],dx[1],dx[2]),
                                          AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                                          parm);

                        for (int v = 0; v < nvar; ++v)
                            accum[v] += wx * Conservative[v];
                    }
#endif

                    for (int v = 0; v < nvar; ++v)
                        dof(i,j,k,v) = accum[v] * inv_norm;
                }

            } // i
        } // j
    } // k
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(amrex::Box const& bx, amrex::Array4<amrex::Real> const& dof,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx_FinestLevel,
         int const& nvar,
         Parm const& parm,
         bool use_quadrature = true,
         int ID_IB = 0)
{
    set_dof_ca(bx, dof, prob_lo, dx, dx_FinestLevel, nvar, parm, use_quadrature, ID_IB);
}







#endif
