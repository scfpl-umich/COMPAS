// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>


// =======================================
// Initial condition
//
// Heat conduction across a stationary contact at x = prob.x0, at the uniform pressure prob.P0:
// the temperature is prob.T_L on the left and prob.T_R on the right, in both phases, and each
// phase density follows from its EOS, rho_k = (A_k P0 + A_k B_k/(1 + A_k))/(cv_k T). With
// prob.layout = 0 phase 1 fills the left and phase 2 the right, with a trace prob.eps of the
// other and a tanh profile of width eta = PhaseField.Eta_Multiplier dx (finest); with
// prob.layout = 1 the mixture is uniform, alpha_1 = prob.alpha_1.
//
// For one gas and a small jump, conduction at constant pressure gives
//   T(x,t) = (T_L + T_R)/2 - (T_L - T_R)/2 erf((x - x0)/(2 sqrt(chi t))),  chi = kappa/(rho cp),
// with cp = gamma cv and rho the density at (P0, (T_L + T_R)/2).
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real ContactDensity(int ID_Phase, amrex::Real T, Parm const& parm)
{
    const amrex::Real A  = (ID_Phase == 1) ? parm.A1   : parm.A2;
    const amrex::Real B  = (ID_Phase == 1) ? parm.B1   : parm.B2;
    const amrex::Real cv = (ID_Phase == 1) ? parm.cv_1 : parm.cv_2;
    return (A*parm.P0 + A*B/(1.0 + A))/(cv*T);
}

// The temperature of the linearized constant-pressure solution, from the properties of phase 1
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real ContactTemperature(amrex::Real x, amrex::Real t, Parm const& parm)
{
    const amrex::Real T_m = 0.5*(parm.T_L + parm.T_R);
    const amrex::Real dT  = parm.T_L - parm.T_R;
    if (t <= 0.0) { return (x < parm.x0) ? parm.T_L : parm.T_R; }
    const amrex::Real cp  = parm.cv_1*(1.0 + parm.A1)/parm.A1;
    const amrex::Real chi = parm.kappa_1/(ContactDensity(1, T_m, parm)*cp);
    return T_m - 0.5*dT*std::erf((x - parm.x0)/(2.0*std::sqrt(chi*t)));
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

    amrex::Real Alpha1 = parm.alpha_1;
    if (parm.layout == 0){
        const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier*dx_FinestLevel;
        const amrex::Real H   = (eta > 0.0) ? 0.5*(1.0 - std::tanh((x - parm.x0)/(std::sqrt(2.0)*eta)))
                                            : ((x < parm.x0) ? 1.0 : 0.0);
        Alpha1 = parm.eps + (1.0 - 2.0*parm.eps)*H;
    }
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*ContactDensity(1, T, parm), Alpha2*ContactDensity(2, T, parm),
       AMREX_D_DECL(u, v, w),
                    Alpha1*parm.P0, Alpha2*parm.P0,
                    Alpha1);
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
//   energy              int (alpha_1 rho_1 E_1 + alpha_2 rho_2 E_2) dV
//   V_1, V_2            int alpha_k dV, the phase volumes (heated at constant pressure, a
//                       phase expands)
//   Q_L, Q_R            heat lost by x < x0 and gained by x > x0 per unit cross-section area,
//                       int sum_k alpha_k rho_k cp_k (T_L - T_k) dV / Area and
//                       int sum_k alpha_k rho_k cp_k (T_k - T_R) dV / Area
//   Q_e                 the heat through x0 of the erf solution, kappa |T_L - T_R| sqrt(t/(pi chi))
//   errT_L1, errT_max   int |T - T_e| dV / (V |T_L - T_R|) and max |T - T_e|/|T_L - T_R|, with
//                       T the mixture temperature (sum alpha_k rho_k T_k)/rho and T_e the
//                       erf solution
//   errT1_max, errT2_max  max |T_k - T_e|/|T_L - T_R| over the cells with alpha_k > 0.01
//   P_dev               max |P - P0|/P0
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real t = Time(o);
    const auto v = Measure(o, {Sum, Sum, Sum, Sum, Sum, Sum, Max, Max, Max, Max, Sum, Sum, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<13>
    {
        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
        Real T_mix;
        Conservative2Temperature(T_mix, Cons, parm);
        amrex::Array<Real, 2> T_k;
        amrex::Array<bool, 2> Present;
        Conservative2PhaseTemperature(T_k, Present, Cons, parm);
        Real P;
        Conservative2Pressure(P, Cons, parm);

        const Real a1   = Cons[INDEX_VolumeFraction1];
        const Real dT   = std::abs(parm.T_L - parm.T_R);
        const Real T_e  = ContactTemperature(c.x[0], t, parm);
        const Real e1   = (a1 > 0.01)       ? std::abs(T_k[0] - T_e)/dT : 0.0;
        const Real e2   = (1.0 - a1 > 0.01) ? std::abs(T_k[1] - T_e)/dT : 0.0;

        // Enthalpy gained (x > x0) or lost (x < x0) since the start, at constant pressure
        const Real cp1  = parm.cv_1*(1.0 + parm.A1)/parm.A1;
        const Real cp2  = parm.cv_2*(1.0 + parm.A2)/parm.A2;
        const bool left = (c.x[0] < parm.x0);
        const Real T0   = left ? parm.T_L : parm.T_R;
        const Real H    = (Present[0] ? Cons[INDEX_Mass1]*cp1*(T_k[0] - T0) : 0.0)
                        + (Present[1] ? Cons[INDEX_Mass2]*cp2*(T_k[1] - T0) : 0.0);
        return {Cons[INDEX_Mass1]*c.dV, Cons[INDEX_Mass2]*c.dV,
                (Cons[INDEX_Energy1] + Cons[INDEX_Energy2])*c.dV,
                a1*c.dV, (1.0 - a1)*c.dV,
                std::abs(T_mix - T_e)/dT*c.dV,
                std::abs(T_mix - T_e)/dT, e1, e2, std::abs(P - parm.P0)/parm.P0,
                c.dV, left ? -H*c.dV/parm.Area : 0.0, left ? 0.0 : H*c.dV/parm.Area};
    });
    Parm const& hp = HostParm(o);
    const Real cp1 = hp.cv_1*(1.0 + hp.A1)/hp.A1;
    const Real chi = hp.kappa_1/(ContactDensity(1, 0.5*(hp.T_L + hp.T_R), hp)*cp1);
    const Real Q_e = hp.kappa_1*std::abs(hp.T_L - hp.T_R)*std::sqrt(t/(PI*chi));
    WriteRow(o, "conduction.csv", {{"mass_1", v[0]}, {"mass_2", v[1]}, {"energy", v[2]},
                                   {"V_1", v[3]}, {"V_2", v[4]},
                                   {"Q_L", v[11]}, {"Q_R", v[12]}, {"Q_e", Q_e}, {"errT_L1", v[5]/v[10]},
                                   {"errT_max", v[6]}, {"errT1_max", v[7]}, {"errT2_max", v[8]},
                                   {"P_dev", v[9]}});
}

#endif  // IC_H_
