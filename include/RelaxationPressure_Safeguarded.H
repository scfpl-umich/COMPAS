// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RelaxationPressure_Safeguarded_H_
#define RelaxationPressure_Safeguarded_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>
#include <cmath>
#include <limits>

// =======================================
// Safeguarded instantaneous pressure relaxation of N phases (six-equation models)
//
// In partial form, every EOS of the six-equation models reads
//   alpha_k rho_k e_k = A_k Pi_k + B_k phi_k + alpha_k rho_k D_k,
// with phi_k = alpha_k - b_k alpha_k rho_k (phi_k = alpha_k without a covolume) and
// Pi_k = phi_k p_k. With the phase energies exchanged by -P d(alpha_k) at the final pressure P
// (P_I implicit; Ndanou et al. 2015, Huang 2026 Eq. 32), the relaxed volume fractions are
//   alpha_k(P) = alpha_k - d_k + c_k/(P + pinf_k),
//   c_k = (Pi_k + phi_k pinf_k)/gamma_k,   d_k = phi_k/gamma_k,
//   gamma_k = (A_k + 1)/A_k,   pinf_k = B_k/(A_k + 1),
// the gamma and P_inf of the stiffened gas and the NASG EOS, and gamma_k = 1 + Gamma_k,
// pinf_k = -p_ref,k/(1 + Gamma_k) for the Mie-Gruneisen EOS. For phi_k > 0, c_k > 0 is
// p_k > -pinf_k, the admissible range of the EOS (a NASG phase past its covolume limit has
// phi_k < 0, and c_k > 0 then means p_k < -pinf_k). c_k and d_k are formed from Pi_k and phi_k,
// so no alpha_k is divided by.
//
// A phase is active when alpha_k > Alpha_Trace (>= 0), alpha_k rho_k > 0 and c_k > 0. The others
// are passive: absent, a trace below Alpha_Trace (Physics.pressure_relaxation_trace, default 0),
// a negative volume fraction or a phase without mass (alpha_k rho_k <= 0) left by the hyperbolic
// step or by the rounding of alpha_2 = 1 - alpha_1 (a massless Mie-Gruneisen trace with zero
// energy sits at p = p_ref, above a positive minimum pressure, and would otherwise take most of
// the cell), or a phase energy below the EOS minimum, for which the relaxation function is not
// defined. The relaxation does not change the masses. A passive phase keeps its alpha_k exactly
// (no volume fraction is clipped or rescaled here: bounds on alpha come from the transport) and
// takes no part in the volume exchange; the caller then resets its energy at the common pressure,
// so it reaches that pressure at constant volume. The active phases share the rest of the volume,
// V = 1 - (sum of the passive alpha_k), so that P solves
//   g(P) = sum_active (c_k/(P + pinf_k) - d_k) + S - V = 0,   S = sum_active alpha_k.
// g is convex and decreasing for P > -pinf_m, pinf_m the smallest pinf_k of the active phases,
// and tends to +infinity at that asymptote, so with two or more active phases it has exactly
// one root there (Huang 2026, Theorems 3.1 and 3.2) when g(+infinity) = S - V - sum d_k < 0,
// which holds in the usual case S = V. At the root each active alpha_k(P) lies above
// alpha_k - d_k > 0, and they add up to V to the round-off of g.
//
// The root is sought in x = P + pinf_m > 0. A bracket search by factors of 4 from a starting
// guess (the phi-weighted mean pressure of the active phases) puts x within a factor 4 to the
// left of the root (g(x) > 0 >= g(4x)); Relax_MaxBracket steps cover the range of the doubles,
// so that the root of a tiny trace next to its asymptote (x down to 1e-300) is found too.
// Newton iterations from there approach the root monotonically from the left, g being convex
// and decreasing; an iterate that leaves the bracket is replaced by a bisection, geometric when
// the bracket spans more than a factor 16. g' is formed as (c_k/y)/y, which does not underflow
// where y^2 would. The iteration stops when the step is below Relax_RelTol x, when g is zero to
// its round-off (Relax_ResTol times the sum of the magnitudes of its terms; the step test alone
// creeps when g' is huge), at an exact root, or after Relax_MaxIter iterations. The active
// volume fractions are alpha_k(P) at the root, as they are.
//
// Return value (the path), with Iter the number of Newton or bisection iterations:
//   Relax_PathSolve     two or more active phases, root found
//   Relax_PathOneActive one active phase, which takes the volume V
//   Relax_PathNoActive  no active phase: every volume fraction is kept
//   Relax_PathNoRoot    g has no root (S - V >= sum d_k, a badly unbalanced input): every volume
//                       fraction is kept
// Relax_PathModel (0) is not returned here; the models use it for their own solution (the
// closed form or the original Newton iteration), kept wherever it is admissible.
// =======================================

constexpr amrex::Real Relax_RelTol     = 1.0e-14;
constexpr amrex::Real Relax_ResTol     = 4.0e-15;
constexpr int         Relax_MaxIter    = 100;
constexpr int         Relax_MaxBracket = 1100;  // 4^1100 spans the doubles, 1e-323 to 1e308

constexpr int Relax_PathModel     = 0;
constexpr int Relax_PathSolve     = 1;
constexpr int Relax_PathOneActive = 2;
constexpr int Relax_PathNoActive  = 3;
constexpr int Relax_PathNoRoot    = 4;

// True for a finite number (not NaN, not infinite)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool Relax_Finite (amrex::Real x) noexcept
{
    return (x == x) && (std::abs(x) <= std::numeric_limits<amrex::Real>::max());
}

// The relaxation coefficients c_k, d_k and pinf_k of one phase, and whether it is active
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool Relax_PhaseCoefficients (amrex::Real& c, amrex::Real& d, amrex::Real& pinf,
                              amrex::Real Alpha, amrex::Real Alpha_rho,
                              amrex::Real Phi, amrex::Real Pi,
                              amrex::Real A, amrex::Real B,
                              amrex::Real Alpha_Trace) noexcept
{
    const amrex::Real gamma = (A + 1.0)/A;
    pinf = B/(A + 1.0);
    c = (Pi + Phi*pinf)/gamma;
    d = Phi/gamma;
    return (Alpha > 0.0) && (Alpha > Alpha_Trace) && (Alpha_rho > 0.0) && (c > 0.0) && Relax_Finite(c) && Relax_Finite(d);
}

template <int N>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int RelaxPressure_Safeguarded (amrex::Array<amrex::Real,N>& AlphaR,
                               amrex::Real& PR,
                               int& Iter,
                               amrex::Array<amrex::Real,N> const& Alpha,
                               amrex::Array<amrex::Real,N> const& Alpha_rho,
                               amrex::Array<amrex::Real,N> const& Phi,
                               amrex::Array<amrex::Real,N> const& Pi,
                               amrex::Array<amrex::Real,N> const& A,
                               amrex::Array<amrex::Real,N> const& B,
                               amrex::Real Alpha_Trace) noexcept
{
    amrex::Array<amrex::Real,N> c, d, pinf, shift;
    amrex::Array<bool,N> active;
    int n_active = 0;
    amrex::Real pinf_m = std::numeric_limits<amrex::Real>::max();
    for (int k = 0; k < N; ++k){
        active[k] = Relax_PhaseCoefficients(c[k], d[k], pinf[k], Alpha[k], Alpha_rho[k], Phi[k], Pi[k],
                                            A[k], B[k], Alpha_Trace);
        if (active[k]){
            ++n_active;
            pinf_m = amrex::min(pinf_m, pinf[k]);
        }
    }
    Iter = 0;
    PR   = 0.0;

    // Every phase starts from its own alpha_k; the passive ones keep it
    amrex::Real S_passive = 0.0;
    for (int k = 0; k < N; ++k){
        AlphaR[k] = Alpha[k];
        if (!active[k]) { S_passive += Alpha[k]; }
    }
    const amrex::Real V = 1.0 - S_passive;

    if (n_active == 0){
        return Relax_PathNoActive;
    }

    if (n_active == 1){
        for (int k = 0; k < N; ++k){
            if (active[k]) { AlphaR[k] = V; }
        }
        return Relax_PathOneActive;
    }

    amrex::Real S = 0.0, D = 0.0, S_Phi = 0.0, S_Pi = 0.0;
    for (int k = 0; k < N; ++k){
        if (active[k]){
            S     += Alpha[k];
            D     += d[k];
            S_Phi += Phi[k];
            S_Pi  += Pi[k];
            shift[k] = pinf[k] - pinf_m;  // >= 0, exactly 0 for the phase that sets the asymptote
        }
    }
    const amrex::Real g_inf = (S - V) - D;
    if (!(g_inf < 0.0)){
        return Relax_PathNoRoot;
    }

    // g(x), the sum of the magnitudes of its terms, and g'(x), x = P + pinf_m
    auto g = [&] (amrex::Real x) noexcept -> amrex::Real
    {
        amrex::Real r = g_inf;
        for (int k = 0; k < N; ++k){
            if (active[k]) { r += c[k]/(x + shift[k]); }
        }
        return r;
    };
    auto g_scale = [&] (amrex::Real x) noexcept -> amrex::Real
    {
        amrex::Real r = std::abs(g_inf);
        for (int k = 0; k < N; ++k){
            if (active[k]) { r += c[k]/(x + shift[k]); }
        }
        return r;
    };
    auto dg = [&] (amrex::Real x) noexcept -> amrex::Real
    {
        amrex::Real r = 0.0;
        for (int k = 0; k < N; ++k){
            if (active[k]) { const amrex::Real y = x + shift[k]; r -= (c[k]/y)/y; }
        }
        return r;
    };

    // Starting guess: the phi-weighted mean pressure of the active phases
    // (x = 1 where it is not positive or not finite; the bracket search covers any start)
    amrex::Real x = (S_Phi > 0.0) ? S_Pi/S_Phi + pinf_m : -1.0;
    if (!(x > 0.0 && Relax_Finite(x))) { x = 1.0; }

    // Bracket: a with g(a) > 0, b = 4a with g(b) <= 0
    amrex::Real gx = g(x);
    amrex::Real a = 0.0;
    amrex::Real b = std::numeric_limits<amrex::Real>::max();
    if (gx > 0.0){
        a = x;
        for (int n = 0; n < Relax_MaxBracket; ++n){
            const amrex::Real x4 = 4.0*x;
            const amrex::Real g4 = g(x4);
            if (!(g4 > 0.0)) { b = x4; break; }
            x = x4; gx = g4; a = x;
        }
    }
    else if (gx < 0.0){
        b = x;
        for (int n = 0; n < Relax_MaxBracket; ++n){
            x *= 0.25;
            gx = g(x);
            if (gx > 0.0) { a = x; break; }
            b = x;
            if (!(gx < 0.0)) { break; }  // exact root or NaN
        }
    }

    const bool bracketed = (gx == 0.0) || ((a > 0.0) && (b < std::numeric_limits<amrex::Real>::max()));
    if (!bracketed){
        // Only for a root beyond the range of the doubles (g_inf < 0 and the asymptote give a sign
        // change); every volume fraction is kept
        return Relax_PathNoRoot;
    }

    if (gx != 0.0 && !(std::abs(gx) <= Relax_ResTol*g_scale(x))){
        // Newton from the left, safeguarded by the bracket (none if the start is a root to round-off)
        for (int it = 0; it < Relax_MaxIter; ++it){
            ++Iter;
            amrex::Real xn = x - gx/dg(x);
            if (!(xn > a && xn < b)){
                xn = (b > 16.0*a) ? std::sqrt(a)*std::sqrt(b) : 0.5*(a + b);
            }
            const amrex::Real step = std::abs(xn - x);
            x  = xn;
            gx = g(x);
            if      (gx > 0.0) { a = x; }
            else if (gx < 0.0) { b = x; }
            else               { break; }
            if (step <= Relax_RelTol*x || std::abs(gx) <= Relax_ResTol*g_scale(x)) { break; }
        }
    }

    PR = x - pinf_m;
    for (int k = 0; k < N; ++k){
        if (active[k]) { AlphaR[k] = Alpha[k] - d[k] + c[k]/(x + shift[k]); }
    }
    return Relax_PathSolve;
}

#endif
