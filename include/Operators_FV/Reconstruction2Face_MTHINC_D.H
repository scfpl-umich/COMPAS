// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Reconstruction2Face_MTHINC_D_H_
#define Reconstruction2Face_MTHINC_D_H_

// Point functions of the multi-dimensional THINC (FiniteVolume.THINC.Multi = 1).

// ln(cosh(x)) without overflow
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
THINC_LnCosh_D(amrex::Real x)
{
    amrex::Real ax = std::abs(x);
    return ax + std::log1p(std::exp(-2.0 * ax)) - std::log(2.0);
}

// Distance-like transform of the volume fraction, in cells (Paula et al. 2023, App. D):
// d = ln(a/(1-a))/(2 beta) + d_corr (2a - 1), d_corr = ln(sinh(beta)/beta)/(2 beta), with a clamped
// to [1e-3, 1 - 1e-3] so that pure cells sit about 1.5 cells from the interface.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
THINC_DistanceTransform_D(amrex::Real a, amrex::Real beta, amrex::Real d_corr)
{
    amrex::Real ac = amrex::min(1.0 - 1.0e-3, amrex::max(1.0e-3, a));
    return std::log(ac / (1.0 - ac)) / (2.0 * beta) + d_corr * (2.0 * ac - 1.0);
}

// Linear multi-dimensional THINC (MTHINC, Ii et al. 2012) on the unit cell X in [-1/2, 1/2]^D:
// H(X) = 0.5 (1 + tanh(beta (n.X + d))), n the unit normal in index space (pointing to larger
// alpha). d is fixed by the cell average with a 3-point tensor Gauss rule, solved by Newton in
// D = tanh(beta d - g) with the shift g of Kumar et al. (2021), which converges monotonically.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
MTHINC_Solve_D(amrex::Real beta,
               amrex::Array<amrex::Real,AMREX_SPACEDIM> const& n,
               amrex::Real Alpha_bar)
{
    const amrex::Real xg[3] = {-0.5*0.7745966692414834, 0.0, 0.5*0.7745966692414834};
    const amrex::Real wg[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};

    constexpr int NG = AMREX_D_TERM(3,*3,*3);
    amrex::Real P[NG];
    amrex::Real W[NG];
    amrex::Real P_min = 1.0e30;
    for (int ig = 0; ig < NG; ig++){
        int gx = ig % 3;
        amrex::Real Pg = n[0] * xg[gx];
        amrex::Real Wg = wg[gx];
#if (AMREX_SPACEDIM > 1)
        int gy = (ig / 3) % 3;
        Pg += n[1] * xg[gy];
        Wg *= wg[gy];
#endif
#if (AMREX_SPACEDIM > 2)
        int gz = ig / 9;
        Pg += n[2] * xg[gz];
        Wg *= wg[gz];
#endif
        P[ig] = beta * Pg;
        W[ig] = Wg;
        P_min = amrex::min(P_min, P[ig]);
    }
    amrex::Real g = -P_min + 1.0e-8;
    for (int ig = 0; ig < NG; ig++){
        P[ig] = std::tanh(P[ig] + g);
    }
    amrex::Real Target = 2.0 * Alpha_bar - 1.0;
    amrex::Real D = -1.0;
    for (int it = 0; it < 40; it++){
        amrex::Real F = -Target;
        amrex::Real dF = 0.0;
        for (int ig = 0; ig < NG; ig++){
            amrex::Real den = 1.0 / (1.0 + P[ig] * D);
            F  += W[ig] * (P[ig] + D) * den;
            dF += W[ig] * (1.0 - P[ig] * P[ig]) * den * den;
        }
        amrex::Real dD = -F / amrex::max(dF, 1.0e-300);
        D = amrex::min(1.0 - 1.0e-15, D + dD);
        if (std::abs(dD) < 1.0e-13){ break; }
    }
    return (std::atanh(D) + g) / beta;
}

// Average of H over the face X_dim = Side/2 (Side = +1 or -1) of the unit cell. The integral along
// the first tangential direction is exact (ln cosh); in 3D the second one uses 3-point Gauss.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
MTHINC_FaceAverage_D(amrex::Real beta,
                     amrex::Array<amrex::Real,AMREX_SPACEDIM> const& n,
                     amrex::Real d,
                     int ID_Dimension,
                     int Side)
{
    amrex::Real Base = 0.5 * amrex::Real(Side) * n[ID_Dimension] + d;
#if (AMREX_SPACEDIM == 1)
    return 0.5 * (1.0 + std::tanh(beta * Base));
#else
    int t1 = (ID_Dimension + 1) % AMREX_SPACEDIM;
    amrex::Real b1 = beta * n[t1];
    auto Average1 = [=] (amrex::Real B) -> amrex::Real
    {
        if (std::abs(b1) < 1.0e-6){
            return 0.5 * (1.0 + std::tanh(beta * B));
        }
        return 0.5 + (THINC_LnCosh_D(beta * B + 0.5 * b1) - THINC_LnCosh_D(beta * B - 0.5 * b1)) / (2.0 * b1);
    };
#if (AMREX_SPACEDIM == 2)
    return Average1(Base);
#else
    int t2 = (ID_Dimension + 2) % AMREX_SPACEDIM;
    const amrex::Real xg[3] = {-0.5*0.7745966692414834, 0.0, 0.5*0.7745966692414834};
    const amrex::Real wg[3] = {5.0/18.0, 8.0/18.0, 5.0/18.0};
    amrex::Real Sum = 0.0;
    for (int ig = 0; ig < 3; ig++){
        Sum += wg[ig] * Average1(Base + n[t2] * xg[ig]);
    }
    return Sum;
#endif
#endif
}

#endif
