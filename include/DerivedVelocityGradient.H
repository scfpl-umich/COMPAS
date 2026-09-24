// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef DERIVED_VELOCITYGRADIENT_H_
#define DERIVED_VELOCITYGRADIENT_H_

/* Derived variables from the velocity gradient: vorticity_mag, q_criterion and lambda2.
 * They are the only derived variables that need neighbor cells, so the callers of
 * compute_derived_var fill one extra ghost cell when one of them is requested. This file
 * uses the model's Conservative2Velocity, so it is included after that header. */

#include <string>

/** True if the derived variable needs the neighbors of the cell, one ghost cell deep */
inline bool DerivedVarNeedsNeighbors (std::string const& name)
{
    return name == "vorticity_mag" || name == "q_criterion" || name == "lambda2";
}

/** Vorticity magnitude, Q criterion and lambda2 criterion at cell (i,j,k), from second-order
 * central differences of the velocity u = rho_u/rho. statein must hold the neighbors of (i,j,k).
 * Q = (|Omega|^2 - |S|^2)/2 and lambda2 is the middle eigenvalue of S^2 + Omega^2, where S and
 * Omega are the symmetric and antisymmetric parts of the velocity gradient. In 2D the velocity
 * gradient has no z components. */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
VelocityGradientInvariants (int i,
                            int j,
                            int k,
                            amrex::Real dx,
                            amrex::Real dy,
                            amrex::Real dz,
                            amrex::Array4<amrex::Real const> const& statein,
                            amrex::Real& vorticity_mag,
                            amrex::Real& q_criterion,
                            amrex::Real& lambda2)
{
    // G[a][b] = d u_a / d x_b
    amrex::Real G[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const amrex::Real h[3] = {dx, dy, dz};

    for (int b = 0; b < AMREX_SPACEDIM; b++){
        const int di = (b == 0), dj = (b == 1), dk = (b == 2);
        amrex::Array<amrex::Real, NSTATE> Um;
        amrex::Array<amrex::Real, NSTATE> Up;
        for (int n = 0; n < NSTATE; n++){
            Um[n] = statein(i-di, j-dj, k-dk, n);
            Up[n] = statein(i+di, j+dj, k+dk, n);
        }
        amrex::Real vm[3] = {0.0, 0.0, 0.0};
        amrex::Real vp[3] = {0.0, 0.0, 0.0};
        Conservative2Velocity(AMREX_D_DECL(vm[0], vm[1], vm[2]), Um);
        Conservative2Velocity(AMREX_D_DECL(vp[0], vp[1], vp[2]), Up);
        for (int a = 0; a < AMREX_SPACEDIM; a++){
            G[a][b] = (vp[a] - vm[a])/(2.0*h[b]);
        }
    }

    const amrex::Real wx = G[2][1] - G[1][2];
    const amrex::Real wy = G[0][2] - G[2][0];
    const amrex::Real wz = G[1][0] - G[0][1];
    vorticity_mag = std::sqrt(wx*wx + wy*wy + wz*wz);

    amrex::Real S[3][3], W[3][3];
    q_criterion = 0.0;
    for (int a = 0; a < 3; a++){
        for (int b = 0; b < 3; b++){
            S[a][b] = 0.5*(G[a][b] + G[b][a]);
            W[a][b] = 0.5*(G[a][b] - G[b][a]);
            q_criterion += 0.5*(W[a][b]*W[a][b] - S[a][b]*S[a][b]);
        }
    }

    // A = S^2 + Omega^2 is symmetric
    amrex::Real A[3][3];
    for (int a = 0; a < 3; a++){
        for (int b = 0; b < 3; b++){
            A[a][b] = 0.0;
            for (int c = 0; c < 3; c++){
                A[a][b] += S[a][c]*S[c][b] + W[a][c]*W[c][b];
            }
        }
    }

    // Eigenvalues of a symmetric 3x3 matrix in closed form (Smith, Comm. ACM 4(4), 1961)
    const amrex::Real p1 = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    amrex::Real e_max, e_min;
    const amrex::Real q = (A[0][0] + A[1][1] + A[2][2])/3.0;
    if (p1 == 0.0){
        e_max = amrex::max(A[0][0], A[1][1], A[2][2]);
        e_min = amrex::min(A[0][0], A[1][1], A[2][2]);
    } else {
        const amrex::Real p2 = (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q)
                             + (A[2][2]-q)*(A[2][2]-q) + 2.0*p1;
        const amrex::Real p = std::sqrt(p2/6.0);
        amrex::Real B[3][3];
        for (int a = 0; a < 3; a++){
            for (int b = 0; b < 3; b++){
                B[a][b] = (A[a][b] - (a == b ? q : 0.0))/p;
            }
        }
        const amrex::Real detB = B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[2][1])
                               - B[0][1]*(B[1][0]*B[2][2] - B[1][2]*B[2][0])
                               + B[0][2]*(B[1][0]*B[2][1] - B[1][1]*B[2][0]);
        const amrex::Real r = amrex::min(amrex::max(0.5*detB, -1.0), 1.0);
        const amrex::Real phi = std::acos(r)/3.0;
        const amrex::Real two_pi_3 = 2.0*3.14159265358979323846/3.0;
        e_max = q + 2.0*p*std::cos(phi);
        e_min = q + 2.0*p*std::cos(phi + two_pi_3);
    }
    lambda2 = 3.0*q - e_max - e_min;
}

#endif
