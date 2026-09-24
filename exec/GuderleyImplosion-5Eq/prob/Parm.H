// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PARM_H_
#define PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuMemory.H>

#include <AMReX_Array.H>
#include <AMReX_MLMG.H>

#include <Physics_Parameter_Compressible5Eq.H>
#include <FiniteVolume_Parameter.H>
#include <PhaseField_Parameter.H>
#include <LinearSystem_Parameter.H>
#include "Guderley.H"

struct Parm
{
    //
    // Physics control parameters
    //
    Physics_Parameter Physics_Parm;

    //
    // Finite Volume parameters
    //
    FiniteVolume_Parameter FiniteVolume_Parm;

    //
    // Phase-Field parameters
    //
    PhaseField_Parameter PhaseField_Parm;

    //
    // Linear solver parameters (Grad_QPF)
    //
    LinearSystem_Parameter LinearSystem_Parm{.max_iter = 10'000, .max_fmg_iter = 10'000, .fixed_iter = 10'000, .verbose = 1};

    //
    // Equation of State (EOS)
    //
    int EOS = 1;

    // Ideal gas (EOS = 0)
    amrex::Real gamma_1;
    amrex::Real gamma_2;

    // Stiffened gas (EOS = 1)
    amrex::Real pinf_1, pinf_2;

    // Mie-Grüneisen (EOS = 2)
    amrex::Real gGamma_1 = 0.0, gGamma_2 = 0.0, pref_1 = 0.0, pref_2 = 0.0, eref_1 = 0.0, eref_2 = 0.0;

    // Thermal properties
    amrex::Real cv_1 = 1.0, cv_2 = 1.0;

    amrex::Real mu_1 = 0.0, mu_2 = 0.0, muB_1 = 0.0, muB_2 = 0.0, kappa_1 = 0.0, kappa_2 = 0.0;

    // Convex EOS fit parameters
    amrex::Real A1, A2, B1, B2, C1, C2, D1, D2;

    //
    // Governing equation parameters
    //
    amrex::Real K = 0.0;
    amrex::Array<amrex::Real, NSTATE> UpperBound;
    amrex::Array<amrex::Real, NSTATE> LowerBound;
    amrex::Array<int,        NSTATE> ID_BoundType;
    amrex::Array<amrex::Real, 2>     p_min;

    // Problem parameters (read from prob.* in the inputs file)
    amrex::Real rho1  = 1.0;   // Density of the outer gas (phase 1)
    amrex::Real rho2  = 5.0;   // Density of the inner gas (phase 2)
    amrex::Real P0    = 0.01;  // Pressure of the gases ahead of the shock
    amrex::Real Mach  = 20.0;  // Initial Mach number of the converging shock
    amrex::Real Rs    = 0.8;   // Initial shock radius
    amrex::Real R0    = 0.5;   // Mean interface radius
    amrex::Real a0    = 0.02;  // Perturbation amplitude
    int         mode  = 8;     // Perturbation mode number

    // Guderley similarity solution behind the shock, computed in DynamicInit()
    Guderley::Table ss;
    amrex::Real tau_c = 0.0;   // time at which the shock would reach the axis

    //
    // Initialization routine
    //
    void Initialize();
    #define DYNAMIC_INIT
    void DynamicInit();
};

AMREX_FORCE_INLINE
void Parm::DynamicInit ()
{
    amrex::ParmParse pp("prob.EOS");

    pp.query("EOS",EOS);
    pp.query("gamma_1",gamma_1);
    pp.query("gamma_2",gamma_2);
    pp.query("pinf_1",pinf_1);
    pp.query("pinf_2",pinf_2);

    pp.query("gGamma_1",gGamma_1);
    pp.query("gGamma_2",gGamma_2);
    pp.query("pref_1",pref_1);
    pp.query("pref_2",pref_2);
    pp.query("eref_1",eref_1);
    pp.query("eref_2",eref_2);

    pp.query("cv_1",cv_1);
    pp.query("cv_2",cv_2);

    pp.query("mu_1",mu_1);
    pp.query("mu_2",mu_2);
    pp.query("muB_1",muB_1);
    pp.query("muB_2",muB_2);
    pp.query("kappa_1",kappa_1);
    pp.query("kappa_2",kappa_2);

    amrex::ParmParse ppp("prob");
    ppp.query("rho1", rho1);
    ppp.query("rho2", rho2);
    ppp.query("P0", P0);
    ppp.query("Mach", Mach);
    ppp.query("Rs", Rs);
    ppp.query("R0", R0);
    ppp.query("a0", a0);
    ppp.query("mode", mode);

    Initialize();

    // Similarity solution for a cylindrical shock (j = 1) in phase 1. The shock speed is
    // Mach times the sound speed ahead of it, and R(t) = Rs (1 - t/tau_c)^alpha, so
    // dR/dt = -alpha Rs/tau_c at t = 0. The table in w reaches the outer boundary at stop_time.
    const amrex::Real alpha = Guderley::Exponent(gamma_1, 1);
    tau_c = alpha * Rs / (Mach * std::sqrt(gamma_1 * P0 / rho1));
    amrex::ParmParse ppg("geometry");
    amrex::Vector<amrex::Real> hi(AMREX_SPACEDIM);
    ppg.getarr("prob_hi", hi, 0, AMREX_SPACEDIM);
    amrex::Real stop_time = 0.0;
    amrex::ParmParse().query("stop_time", stop_time);
    const amrex::Real w_lo = -1.2 * std::pow(Rs, 1.0 / alpha) / tau_c * amrex::max(stop_time - tau_c, 0.0)
                           * std::pow(amrex::min(hi[0], hi[1]), -1.0 / alpha);
    Guderley::Solve(ss, gamma_1, 1, std::exp(2.0), w_lo);
    Guderley::SetScales(ss, Rs, tau_c, rho1);
    amrex::Print() << "Guderley similarity solution: alpha = " << ss.alpha
                   << ", collapse time without the interface tau_c = " << tau_c << "\n";

}

#endif
