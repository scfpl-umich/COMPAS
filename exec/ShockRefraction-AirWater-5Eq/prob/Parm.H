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
#include "PhysicsFunctions.H"

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

    // Problem parameters (read from prob.* in the inputs file), SI units
    amrex::Real rho_air   = 1.17665;    // Air density ahead of the shock (kg/m^3)
    amrex::Real rho_water = 1053.018;   // Water density (kg/m^3)
    amrex::Real P0        = 101325.0;   // Pressure of the air and water at rest (Pa)
    amrex::Real Mach      = 1.46;       // Shock Mach number in air
    amrex::Real beta      = 50.0;       // Angle between the shock and the interface (degrees), 90 - theta_w
    amrex::Real x0        = 0.0;        // Apex of the water wedge on the bottom wall (m)
    amrex::Real x_s       = -2.0e-3;    // Initial position of the shock (m)

    // Post-shock air, from the Rankine-Hugoniot conditions in DynamicInit()
    amrex::Real rho_ps = 0.0;           // density
    amrex::Real u_ps   = 0.0;           // velocity, in +x
    amrex::Real P_ps   = 0.0;           // pressure
    amrex::Real U_s    = 0.0;           // shock speed
    amrex::Real theta_w = 0.0;          // wedge angle between the interface and the wall (radians)

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
    ppp.query("rho_air", rho_air);
    ppp.query("rho_water", rho_water);
    ppp.query("P0", P0);
    ppp.query("Mach", Mach);
    ppp.query("beta", beta);
    ppp.query("x0", x0);
    ppp.query("x_s", x_s);

    Initialize();

    // Rankine-Hugoniot conditions for a shock of Mach number Mach moving in +x into air at rest
    amrex::Real shock[5];
    StiffenedShockConditions(Mach, 0.0, P0, rho_air, gamma_1, pinf_1, shock);
    rho_ps = shock[0];
    u_ps   = shock[1];
    P_ps   = shock[2];
    U_s    = shock[4];
    theta_w = (90.0 - beta) * std::acos(-1.0) / 180.0;

    const amrex::Real c_water = std::sqrt(gamma_2 * (P0 + pinf_2) / rho_water);
    amrex::Print() << "Shock refraction: U_s = " << U_s << " m/s, post-shock air rho = " << rho_ps
                   << " kg/m^3, u = " << u_ps << " m/s, P = " << P_ps << " Pa\n"
                   << "  water sound speed = " << c_water << " m/s, trace speed of the shock on the interface = "
                   << U_s / std::sin(beta * std::acos(-1.0) / 180.0) << " m/s\n";
}

#endif
