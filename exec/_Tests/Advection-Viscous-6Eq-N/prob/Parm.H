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

#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
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
    // Problem parameters (read from prob.* in the inputs file)
    //
    // A drop of phase 2 (radius R, center (xc, yc), density rho_2) in phase 1 (density rho_1),
    // at the uniform pressure P0 and velocity (u0, v0). Any further phase is absent: its
    // volume fraction, partial density and internal energy are exactly 0.
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
