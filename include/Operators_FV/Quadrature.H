// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FACE_QUADRATURE_H_
#define FACE_QUADRATURE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace Quadrature
{

using amrex::Real;

/**
 * Gauss–Legendre quadrature on [-1,1]
 *
 * Provides:
 *   nodes(q)
 *   weights(q)
 *
 */

 /***************************************
 * 2-point Gauss (exact to degree 1)
 ***************************************/
struct GaussLegendre1
{
    static constexpr int npts = 1;

    AMREX_GPU_HOST_DEVICE
    static constexpr Real node (int q)
    {
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE
    static constexpr Real weight (int)
    {
        return 2.0_rt;
    }
};

/***************************************
 * 2-point Gauss (exact to degree 3)
 ***************************************/
struct GaussLegendre2
{
    static constexpr int npts = 2;

    AMREX_GPU_HOST_DEVICE
    static constexpr Real node (int q)
    {
        return (q == 0)
            ? -0.5773502691896257645_rt
            :  0.5773502691896257645_rt;
    }

    AMREX_GPU_HOST_DEVICE
    static constexpr Real weight (int)
    {
        return 1.0_rt;
    }
};


/***************************************
 * 3-point Gauss (exact to degree 5)
 ***************************************/
struct GaussLegendre3
{
    static constexpr int npts = 3;

    AMREX_GPU_HOST_DEVICE
    static constexpr Real node (int q)
    {
        return (q == 0)
            ? -0.7745966692414833770_rt
            : (q == 1)
                ?  0.0_rt
                :  0.7745966692414833770_rt;
    }

    AMREX_GPU_HOST_DEVICE
    static constexpr Real weight (int q)
    {
        return (q == 0 || q == 2)
            ? 0.5555555555555555556_rt   // 5/9
            : 0.8888888888888888889_rt;  // 8/9
    }
};


/***************************************
 * 4-point Gauss (exact to degree 7)
 ***************************************/
struct GaussLegendre4
{
    static constexpr int npts = 4;

    AMREX_GPU_HOST_DEVICE
    static constexpr Real node (int q)
    {
        return (q == 0)
            ? -0.8611363115940525752_rt
            : (q == 1)
                ? -0.3399810435848562648_rt
                : (q == 2)
                    ?  0.3399810435848562648_rt
                    :  0.8611363115940525752_rt;
    }

    AMREX_GPU_HOST_DEVICE
    static constexpr Real weight (int q)
    {
        return (q == 0 || q == 3)
            ? 0.3478548451374538574_rt
            : 0.6521451548625461426_rt;
    }
};

} // namespace Quadrature

#endif
