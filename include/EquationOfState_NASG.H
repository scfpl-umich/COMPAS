// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EquationOfState_NASG_H_
#define EquationOfState_NASG_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Algorithm.H>

// Noble-Abel stiffened gas (prob.EOS.EOS = 3; Le Metayer and Saurel, Phys. Fluids 28, 046102, 2016):
//   rho_k e_k = xi_k (A_k p + B_k) + rho_k q_k,   xi_k = 1 - b_k rho_k,
//   A_k = 1/(gamma_k - 1),  B_k = gamma_k pinf_k/(gamma_k - 1),
// with the covolume b_k in the C slot of the EOS coefficients and q_k in the D slot.
// In partial-density form, alpha_k rho_k e_k = phi_k (A_k p + B_k) + alpha_k rho_k q_k, with
// phi_k = alpha_k - b_k alpha_k rho_k = xi_k alpha_k. The rule of the implementation: every
// alpha_k that multiplies A_k or B_k becomes phi_k.
//
// Covolume bound: v_k > b_k, i.e. xi_k > 0. It is enforced as a floor, xi_k >= NASG_Xi_Min, where
// the covolume factor is formed. The floor is applied only when b_k > 0, and with b_k = 0 the
// helpers below return alpha_k and 1.0 exactly, so the stiffened-gas arithmetic is unchanged.

constexpr amrex::Real NASG_Xi_Min = 1.0e-6;

// Five-equation models: phi = max(alpha - b alpha rho, NASG_Xi_Min alpha); alpha when b = 0
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real NASG_Phi (amrex::Real Alpha, amrex::Real Alpha_rho, amrex::Real b) noexcept
{
    return (b > 0.0) ? amrex::max(Alpha - b*Alpha_rho, NASG_Xi_Min*Alpha) : Alpha;
}

// Six-equation models: xi = max(1 - b rho, NASG_Xi_Min) with rho = alpha rho/alpha; 1.0 when
// b = 0 or alpha <= 0 (an absent phase)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real NASG_Xi (amrex::Real Alpha, amrex::Real Alpha_rho, amrex::Real b) noexcept
{
    return (b > 0.0 && Alpha > 0.0) ? amrex::max(1.0 - b*(Alpha_rho/Alpha), NASG_Xi_Min) : 1.0;
}

// Covolume factor from a phase density (THINC face states): max(1 - b rho, NASG_Xi_Min); 1.0 when b = 0
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real NASG_Xi_rho (amrex::Real rho, amrex::Real b) noexcept
{
    return (b > 0.0) ? amrex::max(1.0 - b*rho, NASG_Xi_Min) : 1.0;
}

#endif
