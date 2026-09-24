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
    // Problem parameters (read from prob.* in the inputs file), SI units
    //
    amrex::Array<amrex::Real, NPHASE> rho;   // densities ahead of the shock
    amrex::Real P0       = 101325.0;         // pressure ahead of the shock
    amrex::Real Mach     = 1.24;             // shock Mach number in material 1
    amrex::Real x_shock  = -0.004;           // initial shock position
    amrex::Real x_12     = 0.0;              // mean position of the interface between materials 1 and 2
    amrex::Real L_2      = 0.010;            // thickness of material 2
    amrex::Real amp      = 0.002;            // perturbation amplitude of both interfaces (0: flat)
    amrex::Real lambda   = 0.060;            // perturbation wavelength
    int         phase_23 = -1;               // -1: second interface in anti-phase with the first, +1: in phase

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

   // Problem parameters. With NPHASE = 4 the fourth material is initially absent
   // (prob/inputs_4materials lists its properties)
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
   ppp.query("x_shock", x_shock);
   ppp.query("x_12", x_12);
   ppp.query("L_2", L_2);
   ppp.query("amp", amp);
   ppp.query("lambda", lambda);
   ppp.query("phase_23", phase_23);

   Initialize();

}

#endif
