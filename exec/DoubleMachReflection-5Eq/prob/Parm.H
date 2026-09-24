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
    amrex::Real rho0    = 1.4;        // Density of the gas at rest ahead of the shock
    amrex::Real P0      = 1.0;        // Pressure ahead of the shock (rho0 = gamma P0 gives c0 = 1)
    amrex::Real Mach    = 10.0;       // Shock Mach number
    amrex::Real theta_w = 30.0;       // Wedge angle in degrees, between the shock normal and the wall
    amrex::Real x0      = 1.0/6.0;    // Wedge tip, where the shock meets the wall at t = 0

    // Shock speed and the state behind the shock, computed in DynamicInit()
    amrex::Real Us   = 0.0;           // Normal speed of the shock
    amrex::Real rhoS = 0.0;           // Density behind the shock
    amrex::Real uS   = 0.0;           // Velocity behind the shock, along the shock normal
    amrex::Real vS   = 0.0;           //   (cos theta_w, -sin theta_w)
    amrex::Real PS   = 0.0;           // Pressure behind the shock
    amrex::Real tan_w = 0.0;          // tan(theta_w)
    amrex::Real cos_w = 1.0;          // cos(theta_w)

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
    ppp.query("rho0", rho0);
    ppp.query("P0", P0);
    ppp.query("Mach", Mach);
    ppp.query("theta_w", theta_w);
    ppp.query("x0", x0);

    Initialize();

    // Normal-shock Rankine-Hugoniot relations for an ideal gas (material 1). The gas behind
    // the shock moves along the shock normal, which is tilted by theta_w below the x axis
    // because the x axis is the wall of the wedge. For Mach 10 and theta_w = 30 degrees this
    // gives (rho, u, v, P) = (8, 7.1447, -4.125, 116.5) as in Woodward and Colella (1984).
    const amrex::Real g  = gamma_1;
    const amrex::Real c0 = std::sqrt(g * P0 / rho0);
    const amrex::Real M2 = Mach * Mach;
    const amrex::Real th = theta_w * std::acos(-1.0) / 180.0;
    Us    = Mach * c0;
    rhoS  = rho0 * (g + 1.0) * M2 / ((g - 1.0) * M2 + 2.0);
    PS    = P0 * (2.0 * g * M2 - (g - 1.0)) / (g + 1.0);
    const amrex::Real uN = Us * (1.0 - rho0 / rhoS);
    uS    =  uN * std::cos(th);
    vS    = -uN * std::sin(th);
    tan_w = std::tan(th);
    cos_w = std::cos(th);
    amrex::Print() << "Double Mach reflection: shock speed " << Us << ", state behind the shock (rho, u, v, P) = ("
                   << rhoS << ", " << uS << ", " << vS << ", " << PS << ")\n";

}

#endif
