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
    // Problem parameters (read from prob.* in the inputs file), nondimensional
    //
    amrex::Array<amrex::Real, NPHASE> rho;   // densities ahead of the shock
    amrex::Real P0      = 1.0;               // pressure ahead of the shock
    amrex::Real Mach    = 2.0;               // Mach number of the converging shock in material 1
    amrex::Real Rs      = 0.92;              // initial shock radius
    amrex::Real R_12    = 0.75;              // mean radius of the interface between materials 1 and 2
    amrex::Real R_23    = 0.5;               // mean radius of the interface between materials 2 and 3
    amrex::Real a_12    = 0.03;              // perturbation amplitude of the outer interface
    amrex::Real a_23    = 0.02;              // perturbation amplitude of the inner interface
    int         mode_12 = 16;                // perturbation mode number of the outer interface (even)
    int         mode_23 = 12;                // perturbation mode number of the inner interface (even)

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

   amrex::ParmParse ppp("prob");
   std::vector<amrex::Real> rhoV;
   ppp.getarr("rho", rhoV);
   if (static_cast<int>(rhoV.size()) < NPHASE || static_cast<int>(gammaV.size()) < NPHASE){
      amrex::Abort("prob.rho and prob.EOS.gamma need one value per material (NPHASE).");
   }
   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
      rho[ID_Phase] = rhoV[ID_Phase];
   }
   ppp.query("P0", P0);
   ppp.query("Mach", Mach);
   ppp.query("Rs", Rs);
   ppp.query("R_12", R_12);
   ppp.query("R_23", R_23);
   ppp.query("a_12", a_12);
   ppp.query("a_23", a_23);
   ppp.query("mode_12", mode_12);
   ppp.query("mode_23", mode_23);
   if (mode_12 % 2 != 0 || mode_23 % 2 != 0){
      amrex::Abort("prob.mode_12 and prob.mode_23 must be even for the quarter plane with symmetry on both axes.");
   }

   Initialize();

}

#endif
