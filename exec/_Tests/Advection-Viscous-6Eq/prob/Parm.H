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
    // Problem parameters (read from prob.* in the inputs file)
    //
    // A drop of phase 2 (radius R, center (xc, yc), density rho_2) in phase 1 (density rho_1),
    // at the uniform pressure P0 and velocity (u0, v0)
    amrex::Real R = 0.4, xc = 0.0, yc = 0.0;
    amrex::Real rho_1 = 0.1, rho_2 = 1.0;
    amrex::Real P0 = 1.0, u0 = 1.0, v0 = 1.0;

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

   amrex::ParmParse pp_prob("prob");
   pp_prob.query("R",R);
   pp_prob.query("xc",xc);
   pp_prob.query("yc",yc);
   pp_prob.query("rho_1",rho_1);
   pp_prob.query("rho_2",rho_2);
   pp_prob.query("P0",P0);
   pp_prob.query("u0",u0);
   pp_prob.query("v0",v0);

   Initialize();

}

#endif
