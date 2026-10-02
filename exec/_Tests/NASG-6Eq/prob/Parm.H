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
    // Problem parameters (read from prob.* in the inputs file), in SI units
    //
    // 0: shock tube along x; 1: a drop of material 1 in material 2 moving at uniform pressure and velocity
    int problem = 0;

    // Shock tube: left (x < x_diaphragm) and right states. The defaults are air (material 2)
    // at 1e9 Pa against liquid water (material 1) at 1e5 Pa.
    amrex::Real rho_L = 50.0,   u_L = 0.0, P_L = 1.0e9;
    amrex::Real rho_R = 1050.0, u_R = 0.0, P_R = 1.0e5;
    amrex::Real x_diaphragm = 0.7;

    // Material that fills each side (1 or 2)
    int mat_L = 2, mat_R = 1;

    // Drop: a disk of radius R_drop centered at the origin, both materials at P_drop and T_drop
    // (the phase densities follow from the EOS), moving with velocity (u_drop, u_drop)
    amrex::Real R_drop = 0.4, P_drop = 1.0e5, T_drop = 300.0, u_drop = 50.0;

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
   pp_prob.query("problem",problem);
   pp_prob.query("R_drop",R_drop);
   pp_prob.query("P_drop",P_drop);
   pp_prob.query("T_drop",T_drop);
   pp_prob.query("u_drop",u_drop);
   pp_prob.query("rho_L",rho_L);
   pp_prob.query("u_L",u_L);
   pp_prob.query("P_L",P_L);
   pp_prob.query("rho_R",rho_R);
   pp_prob.query("u_R",u_R);
   pp_prob.query("P_R",P_R);
   pp_prob.query("x_diaphragm",x_diaphragm);
   pp_prob.query("mat_L",mat_L);
   pp_prob.query("mat_R",mat_R);

   Initialize();

}

#endif
