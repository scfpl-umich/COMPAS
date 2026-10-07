// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RELAXATION_6EQ_N_H_
#define RELAXATION_6EQ_N_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>

#include <array>
#include <fstream>
#include <limits>
#include <vector>


// =======================================
// Relaxation test, N-phase six-equation model (NPHASE = 3), SI units
//
// A row of RelaxTest_NStates cells along x, one test state per cell (cell i holds state i, the
// same in every row in y), at rest and out of pressure equilibrium. Phase 1 is a gas, phase 2 a
// liquid (water), phase 3 a second liquid; their EOS (prob.EOS.*) are set in prob/inputs and
// changed by the variants (gamma from 1.19 to 7.15, pinf up to 1e9 Pa, NASG, Mie-Gruneisen).
// The states are those the pressure relaxation has to handle (RelaxTest_State below): the two
// cells of the water-air tube that made the old iteration fail, traces of each phase in the
// others, compressed and in tension, absent phases, a negative volume fraction, phase energies at
// and below the EOS minimum, two liquids in tension, pressure ratios up to 1e7, a gas near vacuum,
// the NASG covolume limit, and mixture energies off the sum of the phase energies (the relaxation
// resets the phase energies from rho E). States 24 to 31 come from the review of the relaxation:
// a gas at zero energy with a liquid trace (the Newton result lost 71 % of the volume), a 1e-200
// gas trace in a liquid in tension (root next to the gas asymptote, beyond the old bracket search),
// cavitation of a gas trace in a liquid near its EOS minimum, and a trace at a positive EOS minimum
// (Mie-Gruneisen variant). States 32 to 39 test the NASG covolume limits 1/b_k of the liquids
// (1/6.61e-4 kg/m^3 for a liquid without a covolume): a water trace past the limit with almost no
// internal energy, as the SIXEQS air-water tube leaves it next to the initial mixed cell, compressed
// liquid and gas traces next to the limit, a liquid 1e-5 from it, an equilibrium next to it
// (unchanged), both liquids next to their limits, a liquid past it with a gas trace, and a mixture
// energy 5 % off the phase energies next to it.
//
// relax.csv (step 0): each state relaxed once, directly with the relaxation of the solver
// (PressureRelaxation_Cell), one line per state:
//   state, label        the state index and its name
//   path, iter          the path (0 Newton result kept, 1 safeguarded solve, 2 one active
//                       phase, 3 none, 4 no root; RelaxationPressure_Safeguarded.H) and the
//                       number of iterations
//   dmass               max_k |change of alpha_k rho_k|/alpha_k rho_k (0: masses are untouched)
//   dmom, drhoE         |change of rho u| and |change of rho E|/|rho E| (0: untouched)
//   de_sum              |sum_k alpha_k rho_k e_k - (rho E - rho |u|^2/2)|/|rho E - rho |u|^2/2|
//                       after the relaxation (the phase energies add up to the mixture energy)
//   sum_err             |sum_k alpha_k - 1|
//   alpha_min           min_k alpha_k (not checked: the relaxation clips no volume fraction, so a
//                       state out of [0, 1] stays out where the Newton result is kept)
//   alpha_rho_min       min_k alpha_k rho_k
//   dP_after, dP_before max_k |p_k - P|/(|P| + max_k |pinf_k|) over the present phases,
//                       P = sum_k alpha_k p_k, after and before the relaxation (pressure equality)
//   P                   the relaxed pressure
//   de_rel              max_k |E_kR - E_k + P (alpha_kR - alpha_k)|/(sum_k |E_k| + |E_kR| + 2 |P|),
//                       E_k = alpha_k rho_k e_k: the energy exchange -P d(alpha_k) of the relaxation
//                       with the implicit P_I (Huang 2026, Eq. 32), checked on the states it holds for
//                       (exact = 1: every present phase above its EOS minimum, rho E equal to the sum
//                       of the phase energies, no trace threshold). The pressure equality alone holds
//                       for any volume fractions, since the phase energies are reset at one pressure.
//   P_adm               min_k (P + pinf_k)/(|P| + max_k |pinf_k|) over the phases with volume and
//                       mass after the relaxation: the relaxed pressure above every asymptote, which
//                       the relaxation keeps for an admissible state (Huang 2026, Thm 3.1); checked
//                       > -1e-14 where exact = 1 (an unrelaxed trace that takes the common pressure
//                       fails it; the round-off of a relaxed pressure near 0 from E is about 1e-16)
//   xi_min              min_k (1 - b_k alpha_k rho_k/alpha_k) after the relaxation over the phases with
//                       a covolume, volume and mass (1 without): the NASG bound v_k > b_k, checked > 0
//                       where exact = 1 (it holds for the relaxed state of an admissible state)
//   n_passive, passive_dalpha  on the safeguarded paths (path 1 to 4), the number of passive phases
//                       (Relax_PhaseCoefficients: absent, a trace, no mass or an energy below the EOS
//                       minimum) and max |alpha_kR - alpha_k| over them, checked = 0: a passive phase
//                       keeps its volume fraction to the bit (0 and 0 on the Newton path)
// The run aborts if a value is out of its tolerance (RelaxTest_Check). Then it takes max_step
// steps with the same states (all of them in one domain, periodic), and insitu.csv records every
// coarse step the integrals of the masses, of rho E and of the momentum, the smallest volume
// fraction and partial density, max |sum_k alpha_k - 1|, the largest pressure difference
// dP_after and the smallest raw covolume factor (phases with alpha_k > 1e-8) over the domain;
// the run aborts on a NaN there. The steps use the first-order scheme: the first hyperbolic step
// sees the unrelaxed states, adjacent and out of the range of their EOS in places, where the
// reconstruction can fail. The mie-gruneisen variant, with a liquid whose
// minimum pressure several states are below, runs the relaxation test only (max_step = 0).
// =======================================

constexpr int RelaxTest_NStates = 40;

struct RelaxTestState
{
    amrex::Real alpha[3], rho[3], p[3];
    amrex::Real dE;  // relative offset of rho E from the sum of the phase internal energies
};

// State s; p_min is the lower bound p > -pinf_k of each EOS, from Parm
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
RelaxTestState RelaxTest_State (int s, Parm const& parm)
{
    const amrex::Real pmin_2 = parm.p_min[1];
    const amrex::Real pmin_3 = parm.p_min[2];
    // Covolume limits of the liquids (that of the NASG water for a liquid without a covolume)
    const amrex::Real rc_2 = (parm.C[1] > 0.0) ? 1.0/parm.C[1] : 1.0/6.61e-4;
    const amrex::Real rc_3 = (parm.C[2] > 0.0) ? 1.0/parm.C[2] : 1.0/6.61e-4;
    switch (s) {
    //                 alpha_1, alpha_2, alpha_3                   rho_1, rho_2, rho_3              p_1, p_2, p_3                 dE
    // Equilibrium gas-liquid mixture
    case  0: return {{0.5, 0.5, 0.0},                             {1.2, 1000.0, 0.0},               {1.0e5, 1.0e5, 0.0},          0.0};
    // The two cells of the water-air tube (Sod-6Eq-N, water gamma 4.4, level 0, first stage) where
    // the old iteration failed: the gas energy is negative (p_1 < -pinf_1 = 0)
    case  1: return {{0.7509795807, 0.2490204193, 0.0},
                     {38.87919017/0.7509795807, 224.2595564/0.2490204193, 0.0},
                     {-976775.9991/0.7509795807, 184147375.8/0.2490204193, 0.0},                                                  -0.006872};
    // (the water volume fraction 1 - alpha_1 of the tube, 0.003473653183 to the digits printed, so
    // that the volume fractions add up to 1 exactly; the partial densities and pressures are those printed)
    case  2: return {{0.9965263468, 1.0 - 0.9965263468, 0.0},
                     {50.77244352/0.9965263468, 1.949498142/(1.0 - 0.9965263468), 0.0},
                     {-136834.019/0.9965263468, 167843.0881/(1.0 - 0.9965263468), 0.0},                                             0.051926};
    // Liquid trace with a drifted density and an energy below the EOS minimum
    case  3: return {{1.0 - 2.0e-5, 2.0e-5, 0.0},                 {1.2, 1400.0, 0.0},               {1.0e5, -9.0e8, 0.0},         0.0};
    // Shock-compressed liquid trace in gas
    case  4: return {{1.0 - 1.0e-6, 1.0e-6, 0.0},                 {1.2, 1400.0, 0.0},               {1.0e5, 1.0e12, 0.0},         0.0};
    // Gas trace in compressed liquid
    case  5: return {{1.0e-6, 1.0 - 1.0e-6, 0.0},                 {1.2, 1000.0, 0.0},               {1.0e5, 1.0e9, 0.0},          0.0};
    // Gas trace in liquid under tension (the root is close to the gas asymptote P = 0)
    case  6: return {{1.0e-8, 1.0 - 1.0e-8, 0.0},                 {1.2, 1000.0, 0.0},               {1.0e5, -1.0e8, 0.0},         0.0};
    // Pressure ratio 1e7
    case  7: return {{0.5, 0.5, 0.0},                             {0.012, 1000.0, 0.0},             {1.0e3, 1.0e10, 0.0},         0.0};
    // Gas near vacuum against compressed liquid
    case  8: return {{0.5, 0.5, 0.0},                             {1.0e-5, 1000.0, 0.0},            {1.0, 1.0e9, 0.0},            0.0};
    // Three phases at three pressures
    case  9: return {{0.3, 0.4, 0.3},                             {1.2, 1000.0, 1200.0},            {1.0e5, 1.0e9, 1.0e7},        0.0};
    // Compressed trace of one liquid in the other, gas absent
    case 10: return {{0.0, 1.0 - 1.0e-6, 1.0e-6},                 {0.0, 1000.0, 1200.0},            {0.0, 1.0e5, 1.0e9},          0.0};
    // Pure liquid, pure gas (the other phases absent)
    case 11: return {{0.0, 1.0, 0.0},                             {0.0, 1000.0, 0.0},               {0.0, 1.0e9, 0.0},            0.0};
    case 12: return {{1.0, 0.0, 0.0},                             {1.2, 0.0, 0.0},                  {1.0e5, 0.0, 0.0},            0.0};
    // Negative gas volume fraction with a positive mass, as a hyperbolic step can leave it
    case 13: return {{-1.0e-10, 1.0 + 1.0e-10, 0.0},              {-0.01, 1000.0, 0.0},             {1.0e5, 1.0e9, 0.0},          0.0};
    // Liquid exactly at its EOS minimum, p_2 = -pinf_2
    case 14: return {{0.5, 0.5, 0.0},                             {1.2, 1000.0, 0.0},               {1.0e5, pmin_2, 0.0},         0.0};
    // Traces at extreme pressures, one below its EOS minimum
    case 15: return {{1.0e-12, 1.0 - 2.0e-12, 1.0e-12},           {1.2, 1000.0, 1200.0},            {1.0e12, 1.0e8, -2.0e9},      0.0};
    // Mixture energy 20 % above and 1 % below the sum of the phase energies
    case 16: return {{0.5, 0.5, 0.0},                             {1.2, 1000.0, 0.0},               {1.0e5, 1.0e5, 0.0},          0.2};
    case 17: return {{0.5, 0.5, 0.0},                             {12.0, 1000.0, 0.0},              {1.0e6, 1.0e8, 0.0},         -0.01};
    // Liquid near (1497) and past (1550 kg/m^3) the covolume limit 1/b = 1513 kg/m^3 of the NASG water
    case 18: return {{0.5, 0.5, 0.0},                             {1.2, 1497.0, 0.0},               {1.0e5, 1.0e9, 0.0},          0.0};
    case 19: return {{0.5, 0.5, 0.0},                             {1.2, 1550.0, 0.0},               {1.0e5, 1.0e8, 0.0},          0.0};
    // Gas and two liquids in equal volumes
    case 20: return {{1.0/3.0, 1.0/3.0, 1.0/3.0},                 {1.2, 1000.0, 1200.0},            {1.0e5, 1.0e9, 1.0e5},        0.0};
    // Second liquid just above its EOS minimum
    case 21: return {{0.1, 0.0, 0.9},                             {1.2, 0.0, 1200.0},               {1.0e5, 0.0, pmin_3 + 5.0e7}, 0.0};
    // Two liquids in tension, gas absent: the root is negative, below the start of the Newton
    // iteration of the model, which then steps past the asymptote
    case 22: return {{0.0, 0.5, 0.5},                             {0.0, 1000.0, 1200.0},            {0.0, -3.0e8, -1.0e8},        0.0};
    // Gas with two compressed liquid traces
    case 23: return {{1.0 - 2.0e-6, 1.0e-6, 1.0e-6},              {1.2, 1000.0, 1200.0},            {1.0e5, 1.0e9, 1.0e10},       0.0};
    // Gas at zero internal energy (p_1 = 0) with a liquid trace: the Newton iteration ends at a
    // subnormal P, where the 1e-300 guard of f hides the gas term, and alpha_1 (1 - 1/gamma_1)
    case 24: return {{1.0 - 1.0e-12, 1.0e-12, 0.0},               {1.2, 1000.0, 0.0},               {0.0, 1.0e5, 0.0},            0.0};
    // Gas trace of 1e-200 in a liquid in tension: the root lies 1e-194 Pa above the gas asymptote
    case 25: return {{1.0e-200, 1.0 - 1.0e-200, 0.0},             {1.2, 1000.0, 0.0},               {1.0e5, -1.0e8, 0.0},         0.0};
    // Gas near vacuum and a liquid in tension: the relaxed pressure is 9e-3 Pa
    case 26: return {{0.5, 0.5, 0.0},                             {1.0e-6, 1000.0, 0.0},            {1.0e-2, -3.0e8, 0.0},        0.0};
    // Cavitation: a gas trace in a liquid near its EOS minimum
    case 27: return {{1.0e-8, 1.0 - 1.0e-8, 0.0},                 {1.2, 1000.0, 0.0},               {1.0e5, 0.999*pmin_2, 0.0},   0.0};
    // Trace of the second liquid 1 Pa above its EOS minimum in a gas (a positive minimum, 1.4e7 Pa,
    // in the Mie-Gruneisen variant: the root is then next to that asymptote)
    case 28: return {{1.0 - 1.0e-15, 0.0, 1.0e-15},               {1.2, 0.0, 1200.0},               {1.0e5, 0.0, pmin_3 + 1.0},   0.0};
    // Three phases, the gas at zero internal energy
    case 29: return {{0.3, 0.4, 0.3},                             {1.2, 1000.0, 1200.0},            {0.0, 1.0e5, 1.0e5},          0.0};
    // Pure liquid in tension (the other phases absent)
    case 30: return {{0.0, 1.0, 0.0},                             {0.0, 1000.0, 0.0},               {0.0, -2.86e8, 0.0},          0.0};
    // Gas trace in two liquids in tension
    case 31: return {{1.0e-10, 0.5, 0.5 - 1.0e-10},               {1.2, 1000.0, 1200.0},            {1.0e5, -3.0e8, -1.0e8},      0.0};
    // Water trace 3 % past the covolume limit in a compressed gas: the floored EOS gives it almost no
    // internal energy (the state of the SIXEQS air-water tube next to the initial mixed cell)
    case 32: return {{1.0 - 2.3e-5, 2.3e-5, 0.0},                 {50.0, 1.0336*rc_2, 0.0},         {9.9e8, 9.9e8, 0.0},          0.0};
    // Compressed water trace next to the covolume limit (xi = 1e-3)
    case 33: return {{1.0 - 1.0e-6, 1.0e-6, 0.0},                 {1.2, 0.999*rc_2, 0.0},           {1.0e5, 1.0e10, 0.0},         0.0};
    // Gas trace in water compressed next to the covolume limit (xi = 1e-2)
    case 34: return {{1.0e-6, 1.0 - 1.0e-6, 0.0},                 {1.2, 0.99*rc_2, 0.0},            {1.0e5, 2.0e10, 0.0},         0.0};
    // Water 1e-5 from the covolume limit against a gas
    case 35: return {{0.5, 0.5, 0.0},                             {1.2, (1.0 - 1.0e-5)*rc_2, 0.0},  {1.0e5, 1.0e9, 0.0},          0.0};
    // Equilibrium next to the covolume limit (xi = 1e-3): unchanged
    case 36: return {{0.5, 0.5, 0.0},                             {50.0, 0.999*rc_2, 0.0},          {1.0e9, 1.0e9, 0.0},          0.0};
    // Gas and both liquids next to their covolume limits
    case 37: return {{0.2, 0.4, 0.4},                             {1.2, 0.999*rc_2, 0.999*rc_3},    {1.0e5, 1.0e9, 5.0e9},        0.0};
    // Water 2 % past the covolume limit with a gas trace
    case 38: return {{1.0e-6, 1.0 - 1.0e-6, 0.0},                 {1.2, 1.02*rc_2, 0.0},            {1.0e5, 1.0e9, 0.0},          0.0};
    // Water next to the covolume limit (xi = 5e-3), mixture energy 5 % above the phase energies
    default: return {{0.5, 0.5, 0.0},                             {1.2, 0.995*rc_2, 0.0},           {1.0e5, 5.0e9, 0.0},          0.05};
    }
}

// True when the relaxation must satisfy E_kR - E_k = -P (alpha_kR - alpha_k) for every phase: every
// present phase above its EOS minimum (and below the NASG covolume limit), absent phases empty,
// volume fractions in [0, 1], rho E equal to the sum of the phase energies, no trace threshold
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool RelaxTest_Exact (int s, Parm const& parm)
{
    const RelaxTestState st = RelaxTest_State(s, parm);
    bool exact = (st.dE == 0.0) && (parm.Physics_Parm.pressure_relaxation_trace == 0.0);
    for (int k = 0; k < NPHASE; ++k){
        const amrex::Real pinf = parm.B[k]/(parm.A[k] + 1.0);
        if (st.alpha[k] == 0.0) { exact = exact && (st.rho[k] == 0.0) && (st.p[k] == 0.0); }
        else {
            exact = exact && (st.alpha[k] > 0.0) && (st.alpha[k] <= 1.0) && (st.rho[k] > 0.0)
                          && (st.p[k] + pinf > 0.0) && (parm.C[k]*st.rho[k] < 1.0);
        }
    }
    return exact;
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
    amrex::ignore_unused(y, z, AMREX_D_DECL(dx, dy, dz), AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    const int s = amrex::min(amrex::max(static_cast<int>(std::floor(x/dx)), 0), RelaxTest_NStates - 1);
    const RelaxTestState st = RelaxTest_State(s, parm);

    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    amrex::Real P = 0.0;
    for (int k = 0; k < NPHASE; ++k){
        Alpha[k]     = st.alpha[k];
        Alpha_rho[k] = st.alpha[k]*st.rho[k];
        Alpha_P[k]   = st.alpha[k]*st.p[k];
        P += Alpha_P[k];
    }
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(0.0, 0.0, 0.0),
                    P);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
    Conservative[INDEX_Total_Energy] *= (1.0 + st.dE);
}


// Largest |p_k - P| over the phases with alpha_k > 0, P = sum_k alpha_k p_k, relative to the scale
// |P| + max_k |pinf_k| (the round-off of p_k = (alpha_k rho_k e_k - alpha_k B_k)/(A_k alpha_k) is
// relative to |p_k| + |pinf_k|)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real RelaxTest_PressureSpread (amrex::Array<amrex::Real, NSTATE> const& Cons, Parm const& parm)
{
    amrex::Array<amrex::Real, NSTATE> Prim;
    Conservative2Primitive(Prim, Cons, parm);
    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    amrex::Real P;
    AMREX_D_TERM(amrex::Real u;, amrex::Real v;, amrex::Real w;)
    Primitive2State(Alpha, Alpha_rho, Alpha_P, AMREX_D_DECL(u, v, w), P, Prim);
    amrex::ignore_unused(Alpha_rho, AMREX_D_DECL(u, v, w));
    amrex::Real spread = 0.0, scale = std::abs(P);
    for (int k = 0; k < NPHASE; ++k){
        if (Alpha[k] > 0.0){
            spread = amrex::max(spread, std::abs(Alpha_P[k]/Alpha[k] - P));
            scale  = amrex::max(scale, std::abs(P) + std::abs(parm.B[k]/(parm.A[k] + 1.0)));
        }
    }
    return spread/scale;
}


// =======================================
// In-situ diagnostics: relax.csv at step 0, insitu.csv every coarse step (see the top of the file)
// =======================================

constexpr int RelaxTest_NValues = 18;

// Tolerances of the relaxation test; returns a description of the first failure, or ""
inline std::string RelaxTest_Check (std::array<amrex::Real, RelaxTest_NValues> const& v)
{
    if (!(v[2] == 0.0))       { return "the masses changed"; }
    if (!(v[3] == 0.0))       { return "the momentum changed"; }
    if (!(v[4] == 0.0))       { return "rho E changed"; }
    if (!(v[5] <= 1.0e-13))   { return "the phase energies do not add up to rho e"; }
    if (!(v[6] <= 1.0e-14))   { return "sum_k alpha_k differs from 1"; }
    if (!(v[17] == 0.0))      { return "a passive phase changed its volume fraction"; }
    if (!(v[8] >= 0.0))       { return "a partial density is negative"; }
    // The round-off of p_k grows as 1/xi_k near the NASG covolume limit (xi_k = 6e-5 in state 19)
    if (!(v[9] <= 1.0e-10))   { return "the phase pressures differ after the relaxation"; }
    if (!(v[1] <= 100.0))     { return "too many iterations"; }
    if (v[12] == 1.0 && !(v[13] <= 1.0e-12)) { return "the energy exchange is not -P d(alpha_k)"; }
    if (v[12] == 1.0 && !(v[14] > -1.0e-14)) { return "the relaxed pressure is below an EOS minimum"; }
    if (v[12] == 1.0 && !(v[15] > 0.0))      { return "a relaxed phase is past its covolume limit"; }
    return "";
}

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    constexpr int M = RelaxTest_NValues;
    const Real lowest = std::numeric_limits<Real>::lowest();

    if (Step(o) == 0){
        static const char* labels[RelaxTest_NStates] = {
            "equilibrium", "tube diaphragm cell", "tube shock cell", "liquid trace below EOS min",
            "compressed liquid trace", "gas trace in liquid", "gas trace in tension", "pressure ratio 1e7",
            "gas near vacuum", "three phases", "liquid trace in liquid", "pure liquid",
            "pure gas", "negative volume fraction", "liquid at EOS min", "traces at extreme pressures",
            "mixture energy +20%", "mixture energy -1%", "NASG near covolume", "NASG past covolume",
            "three equal volumes", "liquid near EOS min", "two liquids in tension", "two liquid traces in gas",
            "gas at zero energy", "gas trace 1e-200 in tension", "gas near vacuum", "cavitation near EOS min",
            "trace at EOS min + 1 Pa", "three phases with gas at zero energy", "pure liquid in tension",
            "gas trace in two liquids in tension", "water trace past covolume", "compressed water trace near covolume",
            "gas trace in water near covolume", "water 1e-5 from covolume", "equilibrium near covolume",
            "both liquids near covolume", "water past covolume with gas trace", "near covolume with mixture energy +5%"};

        std::vector<std::array<Real, M>> rows(RelaxTest_NStates);
        for (int s = 0; s < RelaxTest_NStates; ++s){
            const auto r = Measure<M>(o, Max, 0,
            [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<M>
            {
                Values<M> v;
                for (int n = 0; n < M; ++n) { v[n] = lowest; }
                if (c.i != s || c.j != 0) { return v; }

                amrex::Array<Real, NSTATE> Cons, ConsR;
                for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
                int Iter, Path;
                PressureRelaxation_Cell(ConsR, Cons, parm, PressureRelaxation_IsNASG(parm), Iter, Path);

                Real dmass = 0.0, sum_alpha = 0.0, alpha_min = 1.0e30, alpha_rho_min = 1.0e30, sum_e = 0.0;
                for (int k = 0; k < NPHASE; ++k){
                    const Real m0 = Cons[INDEX_Mass1 + k], m1 = ConsR[INDEX_Mass1 + k];
                    if (m1 != m0) { dmass = amrex::max(dmass, std::abs(m1 - m0)/amrex::max(std::abs(m0), Real(1.0e-300))); }
                    sum_alpha    += ConsR[INDEX_VolumeFraction1 + k];
                    alpha_min     = amrex::min(alpha_min, ConsR[INDEX_VolumeFraction1 + k]);
                    alpha_rho_min = amrex::min(alpha_rho_min, ConsR[INDEX_Mass1 + k]);
                    sum_e        += ConsR[INDEX_Energy1 + k];
                }
                Real dmom = 0.0, rho = 0.0, rho_u2 = 0.0;
                for (int k = 0; k < NPHASE; ++k) { rho += ConsR[INDEX_Mass1 + k]; }
                for (int d = 0; d < AMREX_SPACEDIM; ++d){
                    dmom   = amrex::max(dmom, std::abs(ConsR[INDEX_MomentumX + d] - Cons[INDEX_MomentumX + d]));
                    rho_u2 += ConsR[INDEX_MomentumX + d]*ConsR[INDEX_MomentumX + d];
                }
                const Real rho_e = ConsR[INDEX_Total_Energy] - 0.5*rho_u2/rho;
                const Real drhoE = std::abs(ConsR[INDEX_Total_Energy] - Cons[INDEX_Total_Energy])/std::abs(Cons[INDEX_Total_Energy]);

                amrex::Array<Real, NSTATE> PrimR;
                Conservative2Primitive(PrimR, ConsR, parm);
                const Real PR = PrimR[NSTATE-1];
                Real de_num = 0.0, de_den = 2.0*std::abs(PR);
                for (int k = 0; k < NPHASE; ++k){
                    const Real E0 = Cons[INDEX_Energy1 + k], E1 = ConsR[INDEX_Energy1 + k];
                    de_num  = amrex::max(de_num, std::abs(E1 - E0 + PR*(ConsR[INDEX_VolumeFraction1 + k] - Cons[INDEX_VolumeFraction1 + k])));
                    de_den += std::abs(E0) + std::abs(E1);
                }
                Real pinf_max = 0.0, p_adm = std::numeric_limits<Real>::max();
                for (int k = 0; k < NPHASE; ++k) { pinf_max = amrex::max(pinf_max, std::abs(parm.B[k]/(parm.A[k] + 1.0))); }
                for (int k = 0; k < NPHASE; ++k){
                    if (ConsR[INDEX_VolumeFraction1 + k] > 0.0 && ConsR[INDEX_Mass1 + k] > 0.0) {
                        p_adm = amrex::min(p_adm, (PR + parm.B[k]/(parm.A[k] + 1.0))/(std::abs(PR) + pinf_max));
                    }
                }
                // Raw covolume factors after the relaxation
                Real xi_min = 1.0;
                for (int k = 0; k < NPHASE; ++k){
                    const Real aR = ConsR[INDEX_VolumeFraction1 + k], mR = ConsR[INDEX_Mass1 + k];
                    if (parm.C[k] > 0.0 && aR > 0.0 && mR > 0.0) { xi_min = amrex::min(xi_min, 1.0 - parm.C[k]*mR/aR); }
                }
                // Passive phases of the safeguarded relaxation, with the phi_k and Pi_k of
                // PressureRelaxation_Cell
                Real n_passive = 0.0, passive_dalpha = 0.0;
                if (Path != Relax_PathModel){
                    amrex::Array<Real, NSTATE> Prim;
                    Conservative2Primitive(Prim, Cons, parm);
                    amrex::Array<Real, NPHASE> a_in, m_in, Alpha_Ps;
                    Real P_in;
                    AMREX_D_TERM(Real u_in;, Real v_in;, Real w_in;)
                    Primitive2State(a_in, m_in, Alpha_Ps, AMREX_D_DECL(u_in, v_in, w_in), P_in, Prim);
                    amrex::ignore_unused(P_in, AMREX_D_DECL(u_in, v_in, w_in));
                    const bool NASG = PressureRelaxation_IsNASG(parm);
                    for (int k = 0; k < NPHASE; ++k){
                        Real Phi = a_in[k], Pi = Alpha_Ps[k];
                        if (NASG){
                            Phi = a_in[k] - parm.C[k]*m_in[k];
                            Pi  = (Cons[INDEX_Energy1 + k] - m_in[k]*parm.D[k] - Phi*parm.B[k])/parm.A[k];
                        }
                        Real cc, dd, pinf;
                        if (!Relax_PhaseCoefficients(cc, dd, pinf, a_in[k], m_in[k], Phi, Pi, parm.A[k], parm.B[k],
                                                     parm.Physics_Parm.pressure_relaxation_trace)){
                            n_passive += 1.0;
                            passive_dalpha = amrex::max(passive_dalpha,
                                                        std::abs(ConsR[INDEX_VolumeFraction1 + k] - Cons[INDEX_VolumeFraction1 + k]));
                        }
                    }
                }
                return {Real(Path), Real(Iter), dmass, dmom, drhoE, std::abs(sum_e - rho_e)/std::abs(rho_e),
                        std::abs(sum_alpha - 1.0), alpha_min, alpha_rho_min,
                        RelaxTest_PressureSpread(ConsR, parm), RelaxTest_PressureSpread(Cons, parm), PR,
                        Real(RelaxTest_Exact(s, parm) ? 1.0 : 0.0), de_num/de_den, p_adm, xi_min,
                        n_passive, passive_dalpha};
            });
            for (int n = 0; n < M; ++n) { rows[s][n] = r[n]; }
        }

        std::string failure;
        int failed_state = -1;
        for (int s = 0; s < RelaxTest_NStates && failure.empty(); ++s){
            failure = RelaxTest_Check(rows[s]);
            if (!failure.empty()) { failed_state = s; }
        }

        if (amrex::ParallelDescriptor::IOProcessor()){
            amrex::UtilCreateDirectory(Directory(o), 0755);
            std::ofstream out(Directory(o) + "relax.csv", std::ios::trunc);
            out.precision(std::numeric_limits<Real>::max_digits10);
            out << "state,label,path,iter,dmass,dmom,drhoE,de_sum,sum_err,alpha_min,alpha_rho_min,dP_after,dP_before,P,exact,de_rel,P_adm,xi_min,n_passive,passive_dalpha\n";
            for (int s = 0; s < RelaxTest_NStates; ++s){
                out << s << "," << labels[s];
                for (int n = 0; n < M; ++n) { out << "," << rows[s][n]; }
                out << "\n";
            }
        }
        if (!failure.empty()){
            amrex::Abort("Relaxation test: state " + std::to_string(failed_state) + " (" + labels[failed_state] + "): " + failure);
        }
        amrex::Print() << "Relaxation test: all " << RelaxTest_NStates << " states within tolerance (relax.csv)\n";
    }

    // Integrals and extrema over the domain, every coarse step
    constexpr int K = NPHASE + 1 + AMREX_SPACEDIM + 5;
    Op ops[K];
    for (int n = 0; n < NPHASE + 1 + AMREX_SPACEDIM; ++n) { ops[n] = Sum; }
    ops[NPHASE + 1 + AMREX_SPACEDIM]     = Min;
    ops[NPHASE + 1 + AMREX_SPACEDIM + 1] = Min;
    ops[NPHASE + 1 + AMREX_SPACEDIM + 2] = Max;
    ops[NPHASE + 1 + AMREX_SPACEDIM + 3] = Max;
    ops[NPHASE + 1 + AMREX_SPACEDIM + 4] = Min;
    const auto r = Measure<K>(o, ops, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<K>
    {
        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Values<K> v;
        Real alpha_min = 1.0e30, alpha_rho_min = 1.0e30, sum_alpha = 0.0;
        for (int k = 0; k < NPHASE; ++k){
            v[k] = Cons[INDEX_Mass1 + k]*c.dV;
            alpha_min     = amrex::min(alpha_min, Cons[INDEX_VolumeFraction1 + k]);
            alpha_rho_min = amrex::min(alpha_rho_min, Cons[INDEX_Mass1 + k]);
            sum_alpha    += Cons[INDEX_VolumeFraction1 + k];
        }
        v[NPHASE] = Cons[INDEX_Total_Energy]*c.dV;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { v[NPHASE + 1 + d] = Cons[INDEX_MomentumX + d]*c.dV; }
        v[NPHASE + 1 + AMREX_SPACEDIM]     = alpha_min;
        v[NPHASE + 1 + AMREX_SPACEDIM + 1] = alpha_rho_min;
        v[NPHASE + 1 + AMREX_SPACEDIM + 2] = std::abs(sum_alpha - 1.0);
        v[NPHASE + 1 + AMREX_SPACEDIM + 3] = RelaxTest_PressureSpread(Cons, parm);
        // Raw covolume factor of the phases with a covolume and alpha_k > 1e-8
        Real xi = 1.0;
        for (int k = 0; k < NPHASE; ++k){
            const Real a = Cons[INDEX_VolumeFraction1 + k];
            if (parm.C[k] > 0.0 && a > 1.0e-8) { xi = amrex::min(xi, 1.0 - parm.C[k]*Cons[INDEX_Mass1 + k]/a); }
        }
        v[NPHASE + 1 + AMREX_SPACEDIM + 4] = xi;
        return v;
    });
    std::vector<std::pair<std::string, Real>> columns;
    for (int k = 0; k < NPHASE; ++k) { columns.push_back({"mass_" + std::to_string(k + 1), r[k]}); }
    columns.push_back({"energy", r[NPHASE]});
    columns.push_back({"momentum_x", r[NPHASE + 1]});
    columns.push_back({"momentum_y", r[NPHASE + 2]});
    columns.push_back({"alpha_min", r[NPHASE + 1 + AMREX_SPACEDIM]});
    columns.push_back({"alpha_rho_min", r[NPHASE + 1 + AMREX_SPACEDIM + 1]});
    columns.push_back({"sum_err", r[NPHASE + 1 + AMREX_SPACEDIM + 2]});
    columns.push_back({"dP", r[NPHASE + 1 + AMREX_SPACEDIM + 3]});
    columns.push_back({"xi_min", r[NPHASE + 1 + AMREX_SPACEDIM + 4]});
    WriteRow(o, "insitu.csv", columns);
    for (auto const& col : columns){
        if (!(col.second == col.second)) { amrex::Abort("Relaxation test: NaN in " + col.first + " at step " + std::to_string(Step(o))); }
    }
}

#endif  // RELAXATION_6EQ_N_H_
