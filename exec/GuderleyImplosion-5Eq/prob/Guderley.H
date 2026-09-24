// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef GUDERLEY_H_
#define GUDERLEY_H_

// =======================================
// Guderley's similarity solution for a converging shock (Guderley 1942)
//
// A strong shock R(t) = A (-t)^alpha converges on the axis (j = 1) or center (j = 2) at
// t = 0, into gas at rest with density rho0 and negligible pressure. Behind it, with
// xi = r/R(t) >= 1,
//     u = (r/t) U(xi),   c^2 = (r/t)^2 C(xi),   rho = rho0 G(xi),
// and with lambda = ln(xi), Delta = (U - alpha)^2 - C,
//     dU/dlambda   = N/Delta,  N = -U (U-1)(U-alpha) + C ((1+j) U - 2 (1-alpha)/gamma)
//     dC/dlambda   = -C [2 (U-1) + (gamma-1)((1+j) U + dU/dlambda)] / (U - alpha)
//     dlnG/dlambda = -[(1+j) U + dU/dlambda] / (U - alpha)
// from the strong-shock state at xi = 1. The exponent alpha is the eigenvalue for which
// the solution crosses the sonic line Delta = 0 where N = 0 too. It is found here by
// bisection (0.835323 for gamma = 1.4, j = 1), and the profile is tabulated once, on the
// host, when the problem parameters are read.
//
// Far behind the shock, and after the collapse time, the flow is continued in
// w = xi^(-1/alpha) = A^(1/alpha) (-t) r^(-1/alpha), which is linear in t at fixed r and
// passes through zero at the collapse. With u = -A^(1/alpha) r^(1-1/alpha) V(w),
// c^2 = A^(2/alpha) r^(2-2/alpha) Z(w) and rho = rho0 G(w) the equations are regular there:
//     Delta = (wV - alpha)^2 - w^2 Z,  M = -V (wV-1)(wV-alpha) + w Z ((1+j) wV - 2 (1-alpha)/gamma)
//     dV/dw   = -[V^2 (wV-alpha)(1-alpha) + alpha Z ((1+j) wV - 2 (1-alpha)/gamma) - wVZ] / Delta
//     dZ/dw   = -Z [2V (1-alpha) - alpha (gamma-1)((1+j) V + M/Delta)] / (wV - alpha)
//     dlnG/dw = alpha ((1+j) V + M/Delta) / (wV - alpha)
// This gives the exact state of the gas flowing in through the outer boundary until the
// reflected shock reaches it.
// =======================================

#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

namespace Guderley {

constexpr int NTAB = 1024;

struct Table
{
    amrex::Real alpha   = 0.0;
    amrex::Real gamma   = 1.4;
    amrex::Real lam_max = 0.0;              // U, C, lnG cover lambda = ln(xi) in [0, lam_max]
    amrex::Real U[NTAB], C[NTAB], lnG[NTAB];
    amrex::Real w_hi = 0.0, w_lo = 0.0;     // V, Z, lnGw cover w from w_hi = exp(-lam_max/alpha) down to w_lo
    amrex::Real V[NTAB], Z[NTAB], lnGw[NTAB];
    amrex::Real Rs = 0.0, tau_c = 0.0, A1a = 0.0, rho0 = 1.0;   // R(t) = Rs (1 - t/tau_c)^alpha, A1a = A^(1/alpha)
};

namespace detail {

inline void Derivs (const std::array<double,3>& y, double a, double g, int j,
                    std::array<double,3>& dy, double& Delta, double& N)
{
    const double U = y[0], C = y[1];
    Delta = (U - a) * (U - a) - C;
    N     = -U * (U - 1.0) * (U - a) + C * ((1 + j) * U - 2.0 * (1.0 - a) / g);
    const double DU = N / Delta;
    dy[0] = DU;
    dy[1] = -C * (2.0 * (U - 1.0) + (g - 1.0) * ((1 + j) * U + DU)) / (U - a);
    dy[2] = -((1 + j) * U + DU) / (U - a);
}

inline std::array<double,3> RK4 (const std::array<double,3>& y, double h, double a, double g, int j)
{
    std::array<double,3> k1, k2, k3, k4, t;
    double D, N;
    Derivs(y, a, g, j, k1, D, N);
    for (int n = 0; n < 3; ++n) { t[n] = y[n] + 0.5 * h * k1[n]; }
    Derivs(t, a, g, j, k2, D, N);
    for (int n = 0; n < 3; ++n) { t[n] = y[n] + 0.5 * h * k2[n]; }
    Derivs(t, a, g, j, k3, D, N);
    for (int n = 0; n < 3; ++n) { t[n] = y[n] + h * k3[n]; }
    Derivs(t, a, g, j, k4, D, N);
    std::array<double,3> out;
    for (int n = 0; n < 3; ++n) { out[n] = y[n] + h / 6.0 * (k1[n] + 2.0 * k2[n] + 2.0 * k3[n] + k4[n]); }
    return out;
}

inline std::array<double,3> ShockState (double a, double g)
{
    return {2.0 * a / (g + 1.0), 2.0 * g * (g - 1.0) * a * a / ((g + 1.0) * (g + 1.0)), std::log((g + 1.0) / (g - 1.0))};
}

inline void DerivsW (double w, const std::array<double,3>& y, double a, double g, int j, std::array<double,3>& dy)
{
    const double V = y[0], Z = y[1];
    const double Delta = (w * V - a) * (w * V - a) - w * w * Z;
    const double M = -V * (w * V - 1.0) * (w * V - a) + w * Z * ((1 + j) * w * V - 2.0 * (1.0 - a) / g);
    dy[0] = -(V * V * (w * V - a) * (1.0 - a) + a * Z * ((1 + j) * w * V - 2.0 * (1.0 - a) / g) - w * V * Z) / Delta;
    dy[1] = -Z * (2.0 * V * (1.0 - a) - a * (g - 1.0) * ((1 + j) * V + M / Delta)) / (w * V - a);
    dy[2] = a * ((1 + j) * V + M / Delta) / (w * V - a);
}

// Integrate from the shock toward the sonic line. Returns +1 or -1, the sign of N, if Delta
// reaches zero first, 2 if N changes sign first, and 0 if neither happens. As alpha
// increases the outcome goes from -1 to 2 to +1, and the eigenvalue is at the 2 -> +1 switch.
inline int Approach (double a, double g, int j, double tol, std::vector<std::array<double,4>>* hist)
{
    std::array<double,3> y = ShockState(a, g), dy, dyn;
    double lam = 0.0, h0 = 1.0e-3;
    if (hist) { hist->push_back({lam, y[0], y[1], y[2]}); }
    while (lam < 3.0) {
        double Delta, N;
        Derivs(y, a, g, j, dy, Delta, N);
        if (std::abs(Delta) < tol) { return (N > 0.0) ? 1 : -1; }
        const double dDelta = 2.0 * (y[0] - a) * dy[0] - dy[1];
        double h = h0;
        if (dDelta > 0.0) { h = std::min(h0, 0.5 * std::abs(Delta) / dDelta); }
        std::array<double,3> yn = RK4(y, h, a, g, j);
        double Dn, Nn;
        Derivs(yn, a, g, j, dyn, Dn, Nn);
        if (Dn > 0.0) { h0 *= 0.5; continue; }
        if (N != 0.0 && (Nn > 0.0) != (N > 0.0)) { return 2; }
        y = yn; lam += h;
        if (hist) { hist->push_back({lam, y[0], y[1], y[2]}); }
    }
    return 0;
}

inline double FindAlpha (double g, int j)
{
    // Scan for the bracket, then bisect
    double lo = 0.0, hi = 0.0;
    int prev = 0;
    for (double a = 0.30; a < 0.999; a += 0.005) {
        const int res = Approach(a, g, j, 1.0e-10, nullptr);
        if (res == 1 && prev == 2) { lo = a - 0.005; hi = a; break; }
        prev = res;
    }
    if (hi == 0.0) { amrex::Abort("Guderley: no similarity exponent found for this gamma"); }
    for (int it = 0; it < 55; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (Approach(mid, g, j, 1.0e-10, nullptr) == 1) { hi = mid; } else { lo = mid; }
    }
    return 0.5 * (lo + hi);
}

} // namespace detail

// The similarity exponent alpha for gamma and j
inline double Exponent (double gamma, int j) { return detail::FindAlpha(gamma, j); }

// Solve for alpha, tabulate U, C and ln G on [1, xi_max], and continue in w down to w_lo
inline void Solve (Table& T, double gamma, int j, double xi_max, double w_lo)
{
    using namespace detail;
    const double a = FindAlpha(gamma, j);

    // Up to the sonic point
    std::vector<std::array<double,4>> H;
    Approach(a, gamma, j, 1.0e-5, &H);

    // Cross it along the local slope of the solution, then integrate to xi_max
    const auto& p1 = H[H.size() - 2];
    const auto& p2 = H[H.size() - 1];
    const double dU = p2[1] - p1[1];
    const double step = (dU > 0.0 ? 1.0 : -1.0) * 1.0e-3;
    std::array<double,3> y = {p2[1] + step, p2[2] + (p2[2] - p1[2]) / dU * step, p2[3] + (p2[3] - p1[3]) / dU * step};
    double lam = p2[0] + (p2[0] - p1[0]) / dU * step;
    H.push_back({lam, y[0], y[1], y[2]});
    const double lam_max = std::log(xi_max);
    while (lam < lam_max) {
        y = RK4(y, 1.0e-3, a, gamma, j);
        lam += 1.0e-3;
        H.push_back({lam, y[0], y[1], y[2]});
    }

    // Resample on a uniform lambda grid
    T.alpha = a;
    T.lam_max = lam_max;
    std::size_t m = 0;
    for (int n = 0; n < NTAB; ++n) {
        const double l = lam_max * n / (NTAB - 1);
        while (m + 2 < H.size() && H[m + 1][0] < l) { ++m; }
        const double w = (l - H[m][0]) / (H[m + 1][0] - H[m][0]);
        T.U[n]   = (1.0 - w) * H[m][1] + w * H[m + 1][1];
        T.C[n]   = (1.0 - w) * H[m][2] + w * H[m + 1][2];
        T.lnG[n] = (1.0 - w) * H[m][3] + w * H[m + 1][3];
    }
    T.gamma = gamma;

    // Continue in w from xi_max through the collapse (w = 0) to w_lo
    const double w_hi = std::exp(-lam_max / a);
    T.w_hi = w_hi;
    T.w_lo = std::min(w_lo, 0.5 * w_hi);
    std::array<double,3> z = {T.U[NTAB - 1] / w_hi, T.C[NTAB - 1] / (w_hi * w_hi), T.lnG[NTAB - 1]};
    const int nsub = 16;
    const double h = (T.w_lo - w_hi) / ((NTAB - 1) * nsub);
    double w = w_hi;
    T.V[0] = z[0]; T.Z[0] = z[1]; T.lnGw[0] = z[2];
    for (int n = 1; n < NTAB; ++n) {
        for (int k = 0; k < nsub; ++k) {
            std::array<double,3> k1, k2, k3, k4, t;
            DerivsW(w, z, a, gamma, j, k1);
            for (int q = 0; q < 3; ++q) { t[q] = z[q] + 0.5 * h * k1[q]; }
            DerivsW(w + 0.5 * h, t, a, gamma, j, k2);
            for (int q = 0; q < 3; ++q) { t[q] = z[q] + 0.5 * h * k2[q]; }
            DerivsW(w + 0.5 * h, t, a, gamma, j, k3);
            for (int q = 0; q < 3; ++q) { t[q] = z[q] + h * k3[q]; }
            DerivsW(w + h, t, a, gamma, j, k4);
            for (int q = 0; q < 3; ++q) { z[q] += h / 6.0 * (k1[q] + 2.0 * k2[q] + 2.0 * k3[q] + k4[q]); }
            w += h;
        }
        T.V[n] = z[0]; T.Z[n] = z[1]; T.lnGw[n] = z[2];
    }
}

// Shock position and time scale: R(t) = Rs (1 - t/tau_c)^alpha
inline void SetScales (Table& T, double Rs, double tau_c, double rho0)
{
    T.Rs = Rs; T.tau_c = tau_c; T.rho0 = rho0;
    T.A1a = std::pow(Rs, 1.0 / T.alpha) / tau_c;
}

// U, C and G at xi >= 1
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Evaluate (Table const& T, amrex::Real xi, amrex::Real& U, amrex::Real& C, amrex::Real& G)
{
    const amrex::Real l = std::log(amrex::max(xi, amrex::Real(1.0))) / T.lam_max * (NTAB - 1);
    const int n = amrex::min(static_cast<int>(l), NTAB - 2);
    const amrex::Real w = amrex::min(l - n, amrex::Real(1.0));
    U = (1.0 - w) * T.U[n]   + w * T.U[n + 1];
    C = (1.0 - w) * T.C[n]   + w * T.C[n + 1];
    G = std::exp((1.0 - w) * T.lnG[n] + w * T.lnG[n + 1]);
}

// Write the constants and both tables to similarity.csv, similarity_xi.csv and similarity_w.csv,
// from the I/O rank
inline void Write (Table const& T, double P0, std::string const& dir)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream s(dir + "similarity.csv"), sx(dir + "similarity_xi.csv"), sw(dir + "similarity_w.csv");
    s.precision(12); sx.precision(12); sw.precision(12);
    s << "alpha,tau_c,Rs,A1a,gamma,rho0,P0\n"
      << T.alpha << "," << T.tau_c << "," << T.Rs << "," << T.A1a << "," << T.gamma << "," << T.rho0 << "," << P0 << "\n";
    sx << "xi,U,C,G\n";
    sw << "w,V,Z,G\n";
    for (int n = 0; n < NTAB; ++n) {
        sx << std::exp(T.lam_max * n / (NTAB - 1)) << "," << T.U[n] << "," << T.C[n] << "," << std::exp(T.lnG[n]) << "\n";
        sw << T.w_hi + (T.w_lo - T.w_hi) * n / (NTAB - 1) << "," << T.V[n] << "," << T.Z[n] << "," << std::exp(T.lnGw[n]) << "\n";
    }
}

// Density, radial velocity and pressure of the similarity solution at radius r and time t.
// Returns false ahead of the converging shock, where the gas is at rest.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool State (Table const& T, amrex::Real r, amrex::Real t, amrex::Real& rho, amrex::Real& u_r, amrex::Real& P)
{
    const amrex::Real dt = T.tau_c - t;
    if (dt > 0.0 && r < T.Rs * std::pow(dt / T.tau_c, T.alpha)) { return false; }
    const amrex::Real w = T.A1a * dt * std::pow(r, -1.0 / T.alpha);
    amrex::Real G, c2;
    if (w >= T.w_hi) {
        amrex::Real U, C;
        Evaluate(T, std::pow(w, -T.alpha), U, C, G);
        u_r = -r / dt * U;
        c2  = (r / dt) * (r / dt) * C;
    }
    else {
        const amrex::Real l = (w - T.w_hi) / (T.w_lo - T.w_hi) * (NTAB - 1);
        const int n = amrex::max(0, amrex::min(static_cast<int>(l), NTAB - 2));
        const amrex::Real f = amrex::min(amrex::max(l - n, amrex::Real(0.0)), amrex::Real(1.0));
        const amrex::Real V = (1.0 - f) * T.V[n] + f * T.V[n + 1];
        const amrex::Real Z = (1.0 - f) * T.Z[n] + f * T.Z[n + 1];
        G = std::exp((1.0 - f) * T.lnGw[n] + f * T.lnGw[n + 1]);
        const amrex::Real s = T.A1a * std::pow(r, 1.0 - 1.0 / T.alpha);
        u_r = -s * V;
        c2  = s * s * Z;
    }
    rho = T.rho0 * G;
    P   = rho * c2 / T.gamma;
    return true;
}

} // namespace Guderley

#endif
