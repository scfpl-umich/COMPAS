// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Heat conduction between two materials
// A problem of the Diffusion case (prob.problem = HeatConduction2Mat), included by
// ProblemICBC.H, for the models 5Eq, 6Eq and 6Eq-N (NPHASE = 2).

#ifndef DIFFUSION_HEATCONDUCTION2MAT_H_
#define DIFFUSION_HEATCONDUCTION2MAT_H_

#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
namespace HeatConduction2Mat {

#if (PHYSICS == SIXEQS_IE_NPHASE)
static_assert(NPHASE == 2, "HeatConduction2Mat sets two phases; build it with NPHASE = 2.");
#endif

// =======================================
// Initial condition
//
// Heat conduction between two materials in contact at x = prob.x0, at rest at the uniform
// pressure prob.P0: phase 1 fills the left at the temperature prob.T_L and phase 2 the right
// at prob.T_R, each with a trace prob.eps of the other, with a step (PhaseField.Eta_Multiplier
// = 0) or a tanh profile of width eta = PhaseField.Eta_Multiplier dx (finest). Each phase
// density follows from its EOS, rho_k = (A_k P0 + A_k B_k/(1 + A_k))/(cv_k T).
//
// For a small jump, conduction between two semi-infinite media gives the interface temperature
// T_i = (e_L T_L + e_R T_R)/(e_L + e_R), with the effusivity e = sqrt(kappa rho cp), and
//   T(x,t) = T_i + (T_side - T_i) erf(|x - x0|/(2 sqrt(chi_side t))),  chi = kappa/(rho cp),
// with rho and cp = cv (1 + A)/A of each side at its initial state (Parm.H).
//
// The five-equation model has one temperature per cell, the mixture temperature. The N-phase
// six-equation model has no heat exchange between the phases (pressure relaxation only), and
// a phase conducts only through itself, so no heat crosses the interface: a known limitation
// of the model.
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real ContactDensity(int ID_Phase, amrex::Real T, Parm const& parm)
{
#if (PHYSICS == SIXEQS_IE_NPHASE)
    // ID_Phase from 0
    const amrex::Real A  = parm.A[ID_Phase];
    const amrex::Real B  = parm.B[ID_Phase];
    const amrex::Real cv = parm.cv[ID_Phase];
#else
    // ID_Phase 1 or 2
    const amrex::Real A  = (ID_Phase == 1) ? parm.A1   : parm.A2;
    const amrex::Real B  = (ID_Phase == 1) ? parm.B1   : parm.B2;
    const amrex::Real cv = (ID_Phase == 1) ? parm.cv_1 : parm.cv_2;
#endif
    return (A*parm.P0 + A*B/(1.0 + A))/(cv*T);
}

// The temperature of the exact solution
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real ContactTemperature(amrex::Real x, amrex::Real t, Parm const& parm)
{
    const bool left = (x < parm.x0);
    if (t <= 0.0) { return left ? parm.T_L : parm.T_R; }
    const amrex::Real T_side = left ? parm.T_L   : parm.T_R;
    const amrex::Real chi    = left ? parm.chi_L : parm.chi_R;
    return parm.T_i + (T_side - parm.T_i)*std::erf(std::abs(x - parm.x0)/(2.0*std::sqrt(chi*t)));
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
{
    // Prevent unused-variable warnings
    sharp = false;
    amrex::ignore_unused(y, z, AMREX_D_DECL(dx, dy, dz));

    const amrex::Real T = (x < parm.x0) ? parm.T_L : parm.T_R;

    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier*dx_FinestLevel;
    const amrex::Real H   = (eta > 0.0) ? 0.5*(1.0 - std::tanh((x - parm.x0)/(std::sqrt(2.0)*eta)))
                                        : ((x < parm.x0) ? 1.0 : 0.0);
    const amrex::Real Alpha1 = parm.eps + (1.0 - 2.0*parm.eps)*H;
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

#if (PHYSICS == FIVEEQS)
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*ContactDensity(1, T, parm), Alpha2*ContactDensity(2, T, parm),
       AMREX_D_DECL(u, v, w),
                    parm.P0,
                    Alpha1);
#elif (PHYSICS == SIXEQS)
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*ContactDensity(1, T, parm), Alpha2*ContactDensity(2, T, parm),
       AMREX_D_DECL(u, v, w),
                    Alpha1*parm.P0, Alpha2*parm.P0,
                    Alpha1);
#else
    amrex::Array<amrex::Real, NPHASE> Alpha = {Alpha1, Alpha2};
    amrex::Array<amrex::Real, NPHASE> Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha[i_Phase]*ContactDensity(i_Phase, T, parm);
        Alpha_P  [i_Phase] = Alpha[i_Phase]*parm.P0;
    }
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    parm.P0);
#endif
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/conduction.csv with, over the composite grid,
//   mass_1, mass_2      int alpha_k rho_k dV
//   energy              int rho E dV (int (alpha_1 rho_1 E_1 + alpha_2 rho_2 E_2) dV in SIXEQS)
//   V_1, V_2            int alpha_k dV, the phase volumes
//   Q_L, Q_R            heat lost by x < x0 and gained by x > x0 per unit cross-section area,
//                       from dQ = dH - V dP: with H = sum_k alpha_k rho_k cp_k (T_k - T_side)
//                       - (P - P0) (T_k the mixture temperature in the five-equation model),
//                       Q_L = -int H dV / Area over x < x0 and Q_R = int H dV / Area over
//                       x > x0. The P term removes the work of the uniform pressure change in
//                       the closed box (the two gases change volume differently when heated).
//   Q_e                 the heat through x0 of the exact solution, 2 e_R (T_i - T_R) sqrt(t/pi)
//   T_i                 the exact interface temperature
//   errT_L1, errT_max   over the window |x - x0| < prob.window, which the acoustic waves of the
//                       start do not reach: int |T - T_e| dV / (V_w |T_L - T_R|) and
//                       max |T - T_e|/|T_L - T_R|, with T the mixture temperature of the model
//                       ((sum alpha_k rho_k T_k)/rho in the six-equation models) and T_e the
//                       exact solution
//   P_dev               max |P - P0|/P0
// =======================================

template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real t = Time(o);
    const auto v = Measure(o, {Sum, Sum, Sum, Sum, Sum, Sum, Max, Max, Sum, Sum, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<11>
    {
        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Real T_mix;
        Conservative2Temperature(T_mix, Cons, parm);
#if (PHYSICS == SIXEQS)
        amrex::Array<Real, 2> T_k;
        amrex::Array<bool, 2> Present;
        Conservative2PhaseTemperature(T_k, Present, Cons, parm);
#elif (PHYSICS == SIXEQS_IE_NPHASE)
        amrex::Array<Real, NPHASE> T_k;
        amrex::Array<bool, NPHASE> Present;
        Conservative2PhaseTemperature(T_k, Present, Cons, parm);
#endif
        Real P;
        Conservative2Pressure(P, Cons, parm);

        const Real a1  = Cons[INDEX_VolumeFraction1];
        const Real dT  = std::abs(parm.T_L - parm.T_R);
        const bool win = (std::abs(c.x[0] - parm.x0) < parm.window);
        const Real err = win ? std::abs(T_mix - ContactTemperature(c.x[0], t, parm))/dT : 0.0;

        // Heat gained (x > x0) or lost (x < x0) since the start, dQ = dH - V dP
#if (PHYSICS == SIXEQS_IE_NPHASE)
        const Real cp1  = parm.cv[0]*(1.0 + parm.A[0])/parm.A[0];
        const Real cp2  = parm.cv[1]*(1.0 + parm.A[1])/parm.A[1];
#else
        const Real cp1  = parm.cv_1*(1.0 + parm.A1)/parm.A1;
        const Real cp2  = parm.cv_2*(1.0 + parm.A2)/parm.A2;
#endif
        const bool left = (c.x[0] < parm.x0);
        const Real T0   = left ? parm.T_L : parm.T_R;

        // The masses, energy, second volume fraction and heat H of the cell
#if (PHYSICS == FIVEEQS)
        const Real m1 = Cons[INDEX_Mass1], m2 = Cons[INDEX_Mass2];
        const Real E  = Cons[INDEX_Energy];
        const Real a2 = 1.0 - a1;
        const Real H  = (m1*cp1 + m2*cp2)*(T_mix - T0)
                      - (P - parm.P0);
#elif (PHYSICS == SIXEQS)
        const Real m1 = Cons[INDEX_Mass1], m2 = Cons[INDEX_Mass2];
        const Real E  = Cons[INDEX_Energy1] + Cons[INDEX_Energy2];
        const Real a2 = 1.0 - a1;
        const Real H  = (Present[0] ? m1*cp1*(T_k[0] - T0) : 0.0)
                      + (Present[1] ? m2*cp2*(T_k[1] - T0) : 0.0)
                      - (P - parm.P0);
#else
        const Real m1 = Cons[INDEX_Mass1], m2 = Cons[INDEX_Mass1+1];
        const Real E  = Cons[INDEX_Total_Energy];
        const Real a2 = Cons[INDEX_VolumeFraction1+1];
        const Real H  = (Present[0] ? m1*cp1*(T_k[0] - T0) : 0.0)
                      + (Present[1] ? m2*cp2*(T_k[1] - T0) : 0.0)
                      - (P - parm.P0);
#endif
        return {m1*c.dV, m2*c.dV, E*c.dV, a1*c.dV, a2*c.dV,
                err*c.dV, err, std::abs(P - parm.P0)/parm.P0,
                win ? c.dV : 0.0, left ? -H*c.dV/parm.Area : 0.0, left ? 0.0 : H*c.dV/parm.Area};
    });
    Parm const& hp = HostParm(o);
    const Real Q_e = 2.0*hp.e_R*(hp.T_i - hp.T_R)*std::sqrt(t/PI);
    WriteRow(o, "conduction.csv", {{"mass_1", v[0]}, {"mass_2", v[1]}, {"energy", v[2]},
                                   {"V_1", v[3]}, {"V_2", v[4]},
                                   {"Q_L", v[9]}, {"Q_R", v[10]}, {"Q_e", Q_e}, {"T_i", hp.T_i},
                                   {"errT_L1", v[5]/v[8]}, {"errT_max", v[6]}, {"P_dev", v[7]}});
}

} // namespace HeatConduction2Mat
#endif // DIFFUSION_HAS_HEATCONDUCTION2MAT

#endif // DIFFUSION_HEATCONDUCTION2MAT_H_
