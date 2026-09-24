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

#include <Physics_Parameter_Compressible6Eq.H>

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
    PhaseField_Parameter PhaseField_Parm{.ID_FluxPhase2Flux = 0};

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

    //
    // Problem parameters (read from prob.* in the inputs file), SI units
    //
    amrex::Real r0         = 0.005;       // initial bubble radius
    amrex::Real delta      = -2.0;        // the bubble center is at x = -delta, in units of r0
    amrex::Real rho_liquid = 1000.0;      // liquid density
    amrex::Real rho_gas    = 1.2;         // gas density
    amrex::Real P_liquid   = 10132500.0;  // far-field liquid pressure
    amrex::Real P_gas      = 2339.0;      // gas pressure

    // Units, set in DynamicInit: length r0, time the Rayleigh collapse time t0 (t_c in inputs),
    // velocity u0 = L0/t0, density rho_liquid and pressure P0 = rho0 u0^2
    amrex::Real L0, t0, u0, rho0, P0;

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
   }

   {
      amrex::ParmParse pp("prob");
      pp.query("r0",r0);
      pp.query("delta",delta);
      pp.query("rho_liquid", rho_liquid);
      pp.query("rho_gas", rho_gas);
      pp.query("P_liquid", P_liquid);
      pp.query("P_gas", P_gas);
   }

   // Nondimensionalization, with the Rayleigh collapse time t0 = 0.915 r0 sqrt(rho_liquid/P_liquid)
   rho0 = rho_liquid;
   L0 = r0;
   t0 = 0.915 * r0 * std::sqrt(rho_liquid/P_liquid); 
   u0 = L0/t0;
   P0 = rho0 * u0 * u0;

   // The stiffened-gas pinf is a pressure, so it is scaled too
   pinf_1 /= P0;
   pinf_2 /= P0;

   Initialize();

}

#endif
