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

#include <Physics_Parameter_Compressible5Eq_NPhase.H>
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

    // Convex EOS fit parameters
    amrex::Array<amrex::Real, NPHASE> A, B, C, D;
    amrex::Array<amrex::Real, NPHASE> gamma, pinf, gGamma, pref, eref, cv, mu, muB, kappa;


    //
    // Governing equation parameters
    //
    amrex::Real K = 0.0;
    amrex::Array<amrex::Real, NSTATE> UpperBound;
    amrex::Array<amrex::Real, NSTATE> LowerBound;
    amrex::Array<int,        NSTATE> ID_BoundType;
    amrex::Array<amrex::Real, NPHASE>     p_min;

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

   std::vector<amrex::Real> gammaV(NPHASE), pinfV(NPHASE), gGammaV(NPHASE), prefV(NPHASE), erefV(NPHASE), cvV(NPHASE), muV(NPHASE), muBV(NPHASE), kappaV(NPHASE);

   pp.query("EOS",EOS);
   pp.queryarr("gamma",gammaV);
   pp.queryarr("pinf",pinfV);

   pp.queryarr("gGamma",gGammaV);
   pp.queryarr("pref",prefV);
   pp.queryarr("eref",erefV);

   pp.queryarr("cv",cvV);

   pp.queryarr("mu",muV);
   pp.queryarr("muB",muBV);
   pp.queryarr("kappa",kappaV);

   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
      gamma[ID_Phase] = gammaV[ID_Phase]; 
      pinf[ID_Phase] = pinfV[ID_Phase]; 
      gGamma[ID_Phase] = gGammaV[ID_Phase]; 
      pref[ID_Phase] = prefV[ID_Phase]; 
      eref[ID_Phase] = erefV[ID_Phase]; 
      cv[ID_Phase] = cvV[ID_Phase]; 
      mu[ID_Phase] = muV[ID_Phase]; 
      muB[ID_Phase] = muBV[ID_Phase];
      kappa[ID_Phase] = kappaV[ID_Phase];
   }

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
