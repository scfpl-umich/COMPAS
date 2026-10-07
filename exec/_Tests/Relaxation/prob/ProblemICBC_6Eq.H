// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RELAXATION_6EQ_H_
#define RELAXATION_6EQ_H_

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
// Relaxation test, two-phase six-equation model, SI units
//
// A row of RelaxTest_NStates cells along x, one test state per cell (cell i holds state i, the
// same in every row in y), at rest and out of pressure equilibrium. Material 1 is a gas, material
// 2 a liquid (water); their EOS (prob.EOS.*) are set in prob/inputs and changed by the variants
// (gamma from 1.19 to 6.12, pinf up to 7e8 Pa, NASG, Mie-Gruneisen, pressure-temperature
// relaxation). The states are those the relaxation has to handle (RelaxTest_State below): the two
// cells of the water-air tube that made the N-phase iteration fail, traces of each phase in the
// other, compressed and in tension, absent phases, volume fractions just out of [0, 1], phase
// energies at and below the EOS minimum, pressure ratios up to 1e7, a gas near vacuum and the
// NASG covolume limit. States 24 to 31 come from the review of the relaxation: gas traces in water
// in tension and gases near vacuum, where the relaxed pressure is close to the gas asymptote and
// the root of the quadratic of the closed form loses its digits (volume fractions off by up to
// 0.05), a 1e-200 gas trace (root beyond the old bracket search), a water trace 1 Pa above its EOS
// minimum (a positive minimum in the Mie-Gruneisen variant), and a pure liquid in tension (the old
// 0/0 patch gave the absent gas volume). States 32 to 39 test the NASG covolume limit 1/b_2 of the
// liquid (1/6.61e-4 kg/m^3 for the EOS without a covolume): a liquid trace past the limit with
// almost no internal energy, as the SIXEQS air-water tube leaves it next to the initial mixed cell,
// compressed liquid and gas traces next to the limit, a liquid 1e-5 from it, an equilibrium next to
// it (unchanged), a liquid exactly at it, a liquid past it with a gas trace, and a pressure ratio
// 5e3 next to it.
//
// relax.csv (step 0): each state relaxed once, directly with the relaxation of the solver
// (Relaxation_Instantaneous), one line per state:
//   state, label        the state index and its name
//   path, iter          the path (0 closed form kept, 1 safeguarded solve, 2 one active phase,
//                       3 none, 4 no root, plus 10 where the pressure-temperature relaxation fell
//                       back to the pressure relaxation; RelaxationPressure_Safeguarded.H) and the
//                       number of iterations
//   dmass               max_k |change of alpha_k rho_k|/alpha_k rho_k (0: masses are untouched)
//   dmom                |change of rho u|
//   dE                  |change of E_1 + E_2|/(|E_1| + |E_2| + 1e-4 (|B_1| + |B_2|)) (total energy; the
//                       B_k term measures a state whose energy is a trace's against the round-off of
//                       B_k alpha_k, as alpha_2 = 1 - alpha_1 is stored to an absolute 1e-16)
//   sum_err             |alpha_1 + alpha_2 - 1|
//   alpha_min           min(alpha_1, alpha_2) (not checked: the relaxation clips no volume fraction,
//                       so a state out of [0, 1] stays out where the closed form is kept)
//   alpha_rho_min       min_k alpha_k rho_k
//   dP_after, dP_before max_k |p_k - P|/(|P| + max_k |pinf_k|) over the phases with alpha_k > 1e-8,
//                       P = alpha_1 p_1 + alpha_2 p_2, after and before the relaxation
//   P                   the relaxed pressure
//   de_rel              |E_1R - E_1 + P (alpha_1R - alpha_1)|/(|E_1| + |E_1R| + |E_2| + |E_2R| + 2 |P|):
//                       the energy exchange -P d(alpha_1) of the pressure relaxation with the
//                       implicit P_I (Huang 2026, Eq. 32), checked on the states it holds for (exact
//                       = 1: pressure relaxation, every present phase above its EOS minimum, no trace
//                       threshold). The pressure equality alone holds for any volume fraction, since
//                       the phase energies are reset at one pressure.
//   P_adm               min_k (P + pinf_k)/(|P| + max_k |pinf_k|) over the phases with volume and
//                       mass after the relaxation: the relaxed pressure above every asymptote, which
//                       the relaxation keeps for an admissible state (Huang 2026, Thm 3.1); checked
//                       > -1e-14 where exact = 1 (an unrelaxed trace that takes the common pressure
//                       fails it; the round-off of a relaxed pressure near 0 from E is about 1e-16)
//   xi_min              min_k (1 - b_k alpha_k rho_k/alpha_k) after the relaxation over the phases with
//                       a covolume, volume and mass (1 without): the NASG bound v_k > b_k, checked > 0
//                       where exact = 1 (it holds for the relaxed state of an admissible state)
//   exact_T, dT         max_k |T_k - T_m|/T_m after the relaxation, T_m = max_k T_k, over the phases with
//                       alpha_k > 1e-8 and mass (Conservative2PhaseTemperature), times |P|/(|P| + max_k |pinf_k|):
//                       the temperatures come from the energies reset at P, whose round-off is relative to
//                       |P| + pinf_k (a gas temperature at P << pinf_k has a relative round-off
//                       1e-16 pinf_k/P). Checked <= 1e-12 where exact_T = 1 (pressure-temperature
//                       relaxation, its closed form kept, an admissible state as for exact)
//   n_passive, passive_dalpha  on the safeguarded paths (path 1 to 4, 11 to 14), the number of passive
//                       phases (Relax_PhaseCoefficients: absent, a trace, no mass or an energy below the
//                       EOS minimum) and max |alpha_kR - alpha_k| over them, checked = 0: a passive phase
//                       keeps its volume fraction to the bit (0 and 0 on the other paths)
// The run aborts if a value is out of its tolerance (RelaxTest_Check). Then it takes max_step
// steps with the same states (all of them in one domain, periodic), and insitu.csv records every
// coarse step the integrals of the masses, of E_1 + E_2 and of the momentum, the smallest volume
// fraction and partial density, the largest pressure difference dP_after and the smallest raw
// covolume factor (phases with alpha_k > 1e-8) over the domain; the run aborts on a NaN there.
// The steps use the first-order scheme: the first hyperbolic step sees the unrelaxed states,
// adjacent and out of the range of their EOS in places (volume fractions just out of [0, 1]),
// where the reconstruction can fail. The mie-gruneisen variant, where most
// cells are below the water minimum pressure, runs the relaxation test only (max_step = 0).
// =======================================

constexpr int RelaxTest_NStates = 40;

struct RelaxTestState
{
    amrex::Real alpha_1, rho[2], p[2];
};

// State s; p_min is the lower bound p > -pinf_k of each EOS, from Parm
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
RelaxTestState RelaxTest_State (int s, Parm const& parm)
{
    const amrex::Real pmin_2 = parm.p_min[1];
    // Covolume limit of the liquid (that of the NASG water when b_2 = 0)
    const amrex::Real rc_2 = (parm.C2 > 0.0) ? 1.0/parm.C2 : 1.0/6.61e-4;
    switch (s) {
    //                 alpha_1                    rho_1, rho_2                                   p_1, p_2
    // Equilibrium gas-liquid mixture
    case  0: return {0.5,                         {1.2, 1000.0},                                 {1.0e5, 1.0e5}};
    // The two cells of the water-air tube (Sod-6Eq-N, water gamma 4.4, level 0, first stage) where
    // the N-phase iteration failed: the gas energy is negative (p_1 < -pinf_1 = 0)
    case  1: return {0.7509795807,                {38.87919017/0.7509795807, 224.2595564/0.2490204193},
                                                  {-976775.9991/0.7509795807, 184147375.8/0.2490204193}};
    case  2: return {0.9965263468,                {50.77244352/0.9965263468, 1.949498142/0.003473653183},
                                                  {-136834.019/0.9965263468, 167843.0881/0.003473653183}};
    // Liquid trace with a drifted density and an energy below the EOS minimum
    case  3: return {1.0 - 2.0e-5,                {1.2, 1400.0},                                 {1.0e5, -9.0e8}};
    // Shock-compressed liquid trace in gas
    case  4: return {1.0 - 1.0e-6,                {1.2, 1400.0},                                 {1.0e5, 1.0e12}};
    // Gas trace in compressed liquid
    case  5: return {1.0e-6,                      {1.2, 1000.0},                                 {1.0e5, 1.0e9}};
    // Gas trace in liquid under tension (the root is close to the gas asymptote P = 0)
    case  6: return {1.0e-8,                      {1.2, 1000.0},                                 {1.0e5, -1.0e8}};
    // Pressure ratio 1e7
    case  7: return {0.5,                         {0.012, 1000.0},                               {1.0e3, 1.0e10}};
    // Gas near vacuum against compressed liquid
    case  8: return {0.5,                         {1.0e-5, 1000.0},                              {1.0, 1.0e9}};
    // Liquid trace with a low density (0.15 of nominal)
    case  9: return {1.0 - 1.0e-3,                {1.2, 150.0},                                  {1.0e5, 1.0e5}};
    // Pure liquid, pure gas
    case 10: return {0.0,                         {1.2, 1000.0},                                 {0.0, 1.0e9}};
    case 11: return {1.0,                         {1.2, 1000.0},                                 {1.0e5, 0.0}};
    // Volume fraction just out of [0, 1], with a positive mass, as a hyperbolic step can leave it
    case 12: return {-1.0e-10,                    {-0.01, 1000.0},                               {1.0e5, 1.0e9}};
    case 13: return {1.0 + 1.0e-10,               {1.2, -1000.0},                                {1.0e5, 1.0e9}};
    // Liquid exactly at its EOS minimum, p_2 = -pinf_2
    case 14: return {0.5,                         {1.2, 1000.0},                                 {1.0e5, pmin_2}};
    // Gas trace at an extreme pressure
    case 15: return {1.0e-12,                     {1.2, 1000.0},                                 {1.0e12, 1.0e8}};
    // Liquid near (1497) and past (1550 kg/m^3) the covolume limit 1/b = 1513 kg/m^3 of the NASG water
    case 16: return {0.5,                         {1.2, 1497.0},                                 {1.0e5, 1.0e9}};
    case 17: return {0.5,                         {1.2, 1550.0},                                 {1.0e5, 1.0e8}};
    // Liquid just above its EOS minimum, and in tension
    case 18: return {0.5,                         {1.2, 1000.0},                                 {1.0e5, pmin_2 + 5.0e7}};
    case 19: return {0.5,                         {1.2, 1000.0},                                 {1.0e5, -3.0e8}};
    // Tiny traces
    case 20: return {1.0e-12,                     {1.2, 1000.0},                                 {1.0e5, 1.0e5}};
    case 21: return {1.0 - 1.0e-12,               {1.2, 1000.0},                                 {1.0e5, 1.0e9}};
    // Pressure ratio 1e4 in equal volumes
    case 22: return {0.5,                         {1.2, 1000.0},                                 {1.0e5, 1.0e9}};
    // Hot compressed gas against liquid
    case 23: return {0.9,                         {1.2, 1000.0},                                 {1.0e8, 1.0e5}};
    // Gas trace of 1e-200 in water in tension: the root lies 1e-194 Pa above the gas asymptote
    case 24: return {1.0e-200,                    {1.2, 1000.0},                                 {1.0e5, -1.0e8}};
    // Gas near vacuum and water in tension: the relaxed pressure is 9e-3 Pa
    case 25: return {0.5,                         {1.0e-6, 1000.0},                              {1.0e-2, -3.0e8}};
    // Compressed gas trace in water in tension: the gas expands to 10 % at 2e-8 Pa
    case 26: return {2.8e-19,                     {1.5e-5, 1000.0},                              {7.8e9, -2.75e8}};
    // Cavitation: a gas trace in water near its EOS minimum
    case 27: return {1.0e-8,                      {1.2, 1000.0},                                 {1.0e5, 0.999*pmin_2}};
    // Water trace 1 Pa above its EOS minimum in a gas (a positive minimum, 2.3e7 Pa, in the
    // Mie-Gruneisen variant: the root is then next to that asymptote)
    case 28: return {1.0 - 1.0e-15,               {1.2, 1000.0},                                 {1.0e5, pmin_2 + 1.0}};
    // Pure liquid in tension, the gas absent
    case 29: return {0.0,                         {1.2, 1000.0},                                 {1.0e5, -2.86e8}};
    // Gas at zero internal energy with a water trace
    case 30: return {1.0 - 1.0e-12,               {1.2, 1000.0},                                 {0.0, 1.0e5}};
    // Gas near vacuum at equilibrium, 1 Pa
    case 31: return {0.5,                         {1.0e-5, 1000.0},                              {1.0, 1.0}};
    // Liquid trace 3 % past the covolume limit in a compressed gas: the floored EOS gives it almost no
    // internal energy (the state of the SIXEQS air-water tube next to the initial mixed cell)
    case 32: return {1.0 - 2.3e-5,                {50.0, 1.0336*rc_2},                           {9.9e8, 9.9e8}};
    // Compressed liquid trace next to the covolume limit (xi = 1e-3)
    case 33: return {1.0 - 1.0e-6,                {1.2, 0.999*rc_2},                             {1.0e5, 1.0e10}};
    // Gas trace in a liquid compressed next to the covolume limit (xi = 1e-2)
    case 34: return {1.0e-6,                      {1.2, 0.99*rc_2},                              {1.0e5, 2.0e10}};
    // Liquid 1e-5 from the covolume limit against a gas
    case 35: return {0.5,                         {1.2, (1.0 - 1.0e-5)*rc_2},                    {1.0e5, 1.0e9}};
    // Equilibrium next to the covolume limit (xi = 1e-3): unchanged
    case 36: return {0.5,                         {50.0, 0.999*rc_2},                            {1.0e9, 1.0e9}};
    // Liquid exactly at the covolume limit (xi = 0, floored)
    case 37: return {0.5,                         {1.2, rc_2},                                   {1.0e5, 1.0e9}};
    // Liquid 2 % past the covolume limit with a gas trace
    case 38: return {1.0e-6,                      {1.2, 1.02*rc_2},                              {1.0e5, 1.0e9}};
    // Pressure ratio 5e3 next to the covolume limit (xi = 5e-3)
    default: return {0.5,                         {12.0, 0.995*rc_2},                            {1.0e6, 5.0e9}};
    }
}

// True when the pressure relaxation must satisfy E_1R - E_1 = -P (alpha_1R - alpha_1): pressure
// relaxation, alpha_1 in [0, 1], every present phase above its EOS minimum (and below the NASG
// covolume limit), no trace threshold
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool RelaxTest_Exact (int s, Parm const& parm)
{
    const RelaxTestState st = RelaxTest_State(s, parm);
    const amrex::Real alpha[2] = {st.alpha_1, 1.0 - st.alpha_1};
    const amrex::Real pinf[2]  = {parm.B1/(parm.A1 + 1.0), parm.B2/(parm.A2 + 1.0)};
    const amrex::Real b[2]     = {parm.C1, parm.C2};
    bool exact = !parm.Physics_Parm.pressure_temperature_relaxation
              && (parm.Physics_Parm.pressure_relaxation_trace == 0.0)
              && (st.alpha_1 >= 0.0) && (st.alpha_1 <= 1.0);
    for (int k = 0; k < 2; ++k){
        if (alpha[k] > 0.0){
            exact = exact && (st.rho[k] > 0.0) && (st.p[k] + pinf[k] > 0.0) && (b[k]*st.rho[k] < 1.0);
        }
    }
    return exact;
}

// True when the pressure-temperature relaxation must equalize the phase temperatures: the same
// conditions as RelaxTest_Exact with the pressure-temperature relaxation on (the caller also
// requires its closed form to be kept)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool RelaxTest_ExactT (int s, Parm const& parm)
{
    const RelaxTestState st = RelaxTest_State(s, parm);
    const amrex::Real alpha[2] = {st.alpha_1, 1.0 - st.alpha_1};
    const amrex::Real pinf[2]  = {parm.B1/(parm.A1 + 1.0), parm.B2/(parm.A2 + 1.0)};
    const amrex::Real b[2]     = {parm.C1, parm.C2};
    bool exact = parm.Physics_Parm.pressure_temperature_relaxation
              && (parm.Physics_Parm.pressure_relaxation_trace == 0.0)
              && (st.alpha_1 >= 0.0) && (st.alpha_1 <= 1.0);
    for (int k = 0; k < 2; ++k){
        if (alpha[k] > 0.0){
            exact = exact && (st.rho[k] > 0.0) && (st.p[k] + pinf[k] > 0.0) && (b[k]*st.rho[k] < 1.0);
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
    const amrex::Real Alpha1 = st.alpha_1;
    const amrex::Real Alpha2 = 1.0 - st.alpha_1;

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*st.rho[0], Alpha2*st.rho[1],
       AMREX_D_DECL(0.0, 0.0, 0.0),
                    Alpha1*st.p[0], Alpha2*st.p[1],
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// Largest |p_k - P| over the phases with alpha_k > 1e-8, P = alpha_1 p_1 + alpha_2 p_2, relative to the
// scale |P| + max_k |pinf_k| (the round-off of p_k from the phase energy is relative to |p_k| + |pinf_k|).
// The state holds alpha_1, and alpha_2 = 1 - alpha_1 has an absolute round-off of 1e-16, so the
// pressure of a trace alpha_k < 1e-8 is not defined to better than 1e-16/alpha_k: such traces are skipped.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real RelaxTest_PressureSpread (amrex::Array<amrex::Real, NSTATE> const& Cons, Parm const& parm)
{
    amrex::Array<amrex::Real, NSTATE> Prim;
    Conservative2Primitive(Prim, Cons, parm);
    amrex::Real Alpha1_rho1, Alpha2_rho2, Alpha1_p1, Alpha2_p2, Alpha1;
    AMREX_D_TERM(amrex::Real u;, amrex::Real v;, amrex::Real w;)
    Primitive2State(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u, v, w), Alpha1_p1, Alpha2_p2, Alpha1, Prim);
    amrex::ignore_unused(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u, v, w));
    const amrex::Real Alpha2 = 1.0 - Alpha1;
    const amrex::Real P = Alpha1_p1 + Alpha2_p2;
    amrex::Real spread = 0.0, scale = std::abs(P);
    if (Alpha1 > 1.0e-8){
        spread = amrex::max(spread, std::abs(Alpha1_p1/Alpha1 - P));
        scale  = amrex::max(scale, std::abs(P) + std::abs(parm.B1/(parm.A1 + 1.0)));
    }
    if (Alpha2 > 1.0e-8){
        spread = amrex::max(spread, std::abs(Alpha2_p2/Alpha2 - P));
        scale  = amrex::max(scale, std::abs(P) + std::abs(parm.B2/(parm.A2 + 1.0)));
    }
    return spread/scale;
}


// =======================================
// In-situ diagnostics: relax.csv at step 0, insitu.csv every coarse step (see the top of the file)
// =======================================

constexpr int RelaxTest_NValues = 19;

// Tolerances of the relaxation test; returns a description of the first failure, or ""
inline std::string RelaxTest_Check (std::array<amrex::Real, RelaxTest_NValues> const& v)
{
    if (!(v[2] == 0.0))       { return "the masses changed"; }
    if (!(v[3] == 0.0))       { return "the momentum changed"; }
    // The closed form rebuilds the phase energies from its pressure, so E_1 + E_2 is kept to the
    // round-off of B_k alpha_k, about 1e-13 of |E_1| + |E_2| for a trace of a stiff liquid
    if (!(v[4] <= 1.0e-12))   { return "E_1 + E_2 changed"; }
    if (!(v[5] <= 1.0e-14))   { return "alpha_1 + alpha_2 differs from 1"; }
    if (!(v[18] == 0.0))      { return "a passive phase changed its volume fraction"; }
    if (!(v[7] >= 0.0))       { return "a partial density is negative"; }
    // The round-off of p_k grows as 1/xi_k near the NASG covolume limit
    if (!(v[8] <= 1.0e-10))   { return "the phase pressures differ after the relaxation"; }
    if (!(v[1] <= 100.0))     { return "too many iterations"; }
    if (v[11] == 1.0 && !(v[12] <= 1.0e-12)) { return "the energy exchange is not -P d(alpha_1)"; }
    if (v[11] == 1.0 && !(v[13] > -1.0e-14)) { return "the relaxed pressure is below an EOS minimum"; }
    if (v[11] == 1.0 && !(v[14] > 0.0))      { return "a relaxed phase is past its covolume limit"; }
    if (v[15] == 1.0 && !(v[16] <= 1.0e-12)) { return "the phase temperatures differ after the relaxation"; }
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
            "gas near vacuum", "low-density liquid trace", "pure liquid", "pure gas",
            "negative gas volume fraction", "negative liquid volume fraction", "liquid at EOS min",
            "gas trace at 1e12 Pa", "NASG near covolume", "NASG past covolume", "liquid near EOS min",
            "liquid in tension", "gas trace 1e-12", "liquid trace 1e-12", "pressure ratio 1e4",
            "hot compressed gas", "gas trace 1e-200 in tension", "gas near vacuum", "compressed gas trace in tension",
            "cavitation near EOS min", "trace at EOS min + 1 Pa", "pure liquid in tension",
            "gas at zero energy", "gas at 1 Pa", "liquid trace past covolume", "compressed liquid trace near covolume",
            "gas trace in liquid near covolume", "liquid 1e-5 from covolume", "equilibrium near covolume",
            "liquid at covolume", "liquid past covolume with gas trace", "pressure ratio 5e3 near covolume"};

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
                Relaxation_Instantaneous(ConsR, Cons, parm, Iter, Path);

                Real dmass = 0.0;
                for (int n : {INDEX_Mass1, INDEX_Mass2}){
                    if (ConsR[n] != Cons[n]) { dmass = amrex::max(dmass, std::abs(ConsR[n] - Cons[n])/amrex::max(std::abs(Cons[n]), Real(1.0e-300))); }
                }
                Real dmom = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d){
                    dmom = amrex::max(dmom, std::abs(ConsR[INDEX_MomentumX + d] - Cons[INDEX_MomentumX + d]));
                }
                const Real E  = Cons[INDEX_Energy1] + Cons[INDEX_Energy2];
                const Real ER = ConsR[INDEX_Energy1] + ConsR[INDEX_Energy2];
                const Real a1 = ConsR[INDEX_VolumeFraction1];
                const Real a2 = 1.0 - a1;

                // The relaxed pressure alpha_1 p_1 + alpha_2 p_2 (the primitive state holds alpha_k p_k
                // in the slots of the phase energies)
                amrex::Array<Real, NSTATE> PrimR;
                Conservative2Primitive(PrimR, ConsR, parm);
                const Real PR = PrimR[INDEX_Energy1] + PrimR[INDEX_Energy2];
                const Real de_rel = std::abs(ConsR[INDEX_Energy1] - Cons[INDEX_Energy1] + PR*(a1 - Cons[INDEX_VolumeFraction1]))
                                  /(std::abs(Cons[INDEX_Energy1]) + std::abs(ConsR[INDEX_Energy1]) + std::abs(Cons[INDEX_Energy2])
                                    + std::abs(ConsR[INDEX_Energy2]) + 2.0*std::abs(PR));
                const Real pinf1 = parm.B1/(parm.A1 + 1.0), pinf2 = parm.B2/(parm.A2 + 1.0);
                const Real pinf_max = amrex::max(std::abs(pinf1), std::abs(pinf2));
                Real p_adm = std::numeric_limits<Real>::max();
                if (a1 > 0.0 && ConsR[INDEX_Mass1] > 0.0) { p_adm = amrex::min(p_adm, (PR + pinf1)/(std::abs(PR) + pinf_max)); }
                if (a2 > 0.0 && ConsR[INDEX_Mass2] > 0.0) { p_adm = amrex::min(p_adm, (PR + pinf2)/(std::abs(PR) + pinf_max)); }
                // Raw covolume factors after the relaxation, and the phase temperatures
                Real xi_min = 1.0;
                if (parm.C1 > 0.0 && a1 > 0.0 && ConsR[INDEX_Mass1] > 0.0) { xi_min = amrex::min(xi_min, 1.0 - parm.C1*ConsR[INDEX_Mass1]/a1); }
                if (parm.C2 > 0.0 && a2 > 0.0 && ConsR[INDEX_Mass2] > 0.0) { xi_min = amrex::min(xi_min, 1.0 - parm.C2*ConsR[INDEX_Mass2]/a2); }
                amrex::Array<Real,2> TR;
                amrex::Array<bool,2> PresentR;
                Conservative2PhaseTemperature(TR, PresentR, ConsR, parm);
                const Real aR[2] = {a1, a2};
                Real T_m = 0.0, dT = 0.0;
                for (int k = 0; k < 2; ++k) { if (PresentR[k] && aR[k] > 1.0e-8) { T_m = amrex::max(T_m, std::abs(TR[k])); } }
                for (int k = 0; k < 2; ++k) { if (PresentR[k] && aR[k] > 1.0e-8 && T_m > 0.0) { dT = amrex::max(dT, std::abs(TR[k] - T_m)/T_m); } }
                dT *= std::abs(PR)/(std::abs(PR) + pinf_max);
                const bool exact_T = RelaxTest_ExactT(s, parm) && (Path == Relax_PathModel);
                // Passive phases of the safeguarded relaxation, with the phi_k and Pi_k that the pressure
                // relaxation gives it (the states are at rest)
                Real n_passive = 0.0, passive_dalpha = 0.0;
                if (Path % 10 != Relax_PathModel){
                    amrex::Array<Real, NSTATE> Prim;
                    Conservative2Primitive(Prim, Cons, parm);
                    const Real a_in[2] = {Cons[INDEX_VolumeFraction1], 1.0 - Cons[INDEX_VolumeFraction1]};
                    const Real m_in[2] = {Cons[INDEX_Mass1], Cons[INDEX_Mass2]};
                    const Real A[2] = {parm.A1, parm.A2}, B[2] = {parm.B1, parm.B2};
                    const Real b[2] = {parm.C1, parm.C2}, q[2] = {parm.D1, parm.D2};
                    for (int k = 0; k < 2; ++k){
                        Real Phi = a_in[k], Pi = Prim[INDEX_Energy1 + k];
                        if (parm.C1 > 0.0 || parm.C2 > 0.0){
                            Phi = a_in[k] - b[k]*m_in[k];
                            Pi  = (Cons[INDEX_Energy1 + k] - m_in[k]*q[k] - B[k]*Phi)/A[k];
                        }
                        Real cc, dd, pinf;
                        if (!Relax_PhaseCoefficients(cc, dd, pinf, a_in[k], m_in[k], Phi, Pi, A[k], B[k],
                                                     parm.Physics_Parm.pressure_relaxation_trace)){
                            n_passive += 1.0;
                            passive_dalpha = amrex::max(passive_dalpha, std::abs(aR[k] - a_in[k]));
                        }
                    }
                }
                return {Real(Path), Real(Iter), dmass, dmom,
                        std::abs(ER - E)/(std::abs(Cons[INDEX_Energy1]) + std::abs(Cons[INDEX_Energy2])
                                          + 1.0e-4*(std::abs(parm.B1) + std::abs(parm.B2))),
                        std::abs(a1 + a2 - 1.0), amrex::min(a1, a2),
                        amrex::min(ConsR[INDEX_Mass1], ConsR[INDEX_Mass2]),
                        RelaxTest_PressureSpread(ConsR, parm), RelaxTest_PressureSpread(Cons, parm),
                        PR, Real(RelaxTest_Exact(s, parm) ? 1.0 : 0.0), de_rel, p_adm,
                        xi_min, Real(exact_T ? 1.0 : 0.0), dT, n_passive, passive_dalpha};
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
            out << "state,label,path,iter,dmass,dmom,dE,sum_err,alpha_min,alpha_rho_min,dP_after,dP_before,P,exact,de_rel,P_adm,xi_min,exact_T,dT,n_passive,passive_dalpha\n";
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
    constexpr int K = 3 + AMREX_SPACEDIM + 4;
    Op ops[K] = {Sum, Sum, Sum, AMREX_D_DECL(Sum, Sum, Sum), Min, Min, Max, Min};
    const auto r = Measure<K>(o, ops, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<K>
    {
        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Values<K> v;
        v[0] = Cons[INDEX_Mass1]*c.dV;
        v[1] = Cons[INDEX_Mass2]*c.dV;
        v[2] = (Cons[INDEX_Energy1] + Cons[INDEX_Energy2])*c.dV;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { v[3 + d] = Cons[INDEX_MomentumX + d]*c.dV; }
        v[3 + AMREX_SPACEDIM]     = amrex::min(Cons[INDEX_VolumeFraction1], 1.0 - Cons[INDEX_VolumeFraction1]);
        v[3 + AMREX_SPACEDIM + 1] = amrex::min(Cons[INDEX_Mass1], Cons[INDEX_Mass2]);
        v[3 + AMREX_SPACEDIM + 2] = RelaxTest_PressureSpread(Cons, parm);
        // Raw covolume factor of the phases with a covolume and alpha_k > 1e-8
        const Real a[2] = {Cons[INDEX_VolumeFraction1], 1.0 - Cons[INDEX_VolumeFraction1]};
        const Real b[2] = {parm.C1, parm.C2};
        Real xi = 1.0;
        for (int k = 0; k < 2; ++k){
            if (b[k] > 0.0 && a[k] > 1.0e-8) { xi = amrex::min(xi, 1.0 - b[k]*Cons[INDEX_Mass1 + k]/a[k]); }
        }
        v[3 + AMREX_SPACEDIM + 3] = xi;
        return v;
    });
    std::vector<std::pair<std::string, Real>> columns = {
        {"mass_1", r[0]}, {"mass_2", r[1]}, {"energy", r[2]}, {"momentum_x", r[3]}, {"momentum_y", r[4]},
        {"alpha_min", r[3 + AMREX_SPACEDIM]}, {"alpha_rho_min", r[3 + AMREX_SPACEDIM + 1]},
        {"dP", r[3 + AMREX_SPACEDIM + 2]}, {"xi_min", r[3 + AMREX_SPACEDIM + 3]}};
    WriteRow(o, "insitu.csv", columns);
    for (auto const& col : columns){
        if (!(col.second == col.second)) { amrex::Abort("Relaxation test: NaN in " + col.first + " at step " + std::to_string(Step(o))); }
    }
}

#endif  // RELAXATION_6EQ_H_
