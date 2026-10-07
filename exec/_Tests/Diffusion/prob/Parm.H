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

#include <string>

// One Parm.H for every build of this case (BUILD = 5Eq, 5Eq-N, 6Eq, 6Eq-N, 6Eq-N3 in the
// GNUmakefile) and every problem (prob.problem). The values of FIVEEQS, SIXEQS, ... that PHYSICS
// is compared with are in Macros.H.
#include <Macros.H>
#if (PHYSICS == FIVEEQS)
#include <Physics_Parameter_Compressible5Eq.H>
#elif (PHYSICS == SIXEQS)
#include <Physics_Parameter_Compressible6Eq.H>
#elif (PHYSICS == FIVEEQS_NPHASE)
#include <Physics_Parameter_Compressible5Eq_NPhase.H>
#elif (PHYSICS == SIXEQS_IE_NPHASE)
#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
#endif
#include <FiniteVolume_Parameter.H>
#include <PhaseField_Parameter.H>
#include <LinearSystem_Parameter.H>

// =======================================
// The problems of this case, chosen at run time by prob.problem in the inputs file. A problem
// is compiled only in the builds of the models it is written for; asking for it in another
// build stops the run with a message.
//   Couette2Layer              two-layer Couette flow                   5Eq, 6Eq, 6Eq-N
//   HeatConduction2Mat         heat conduction between two materials    5Eq, 6Eq, 6Eq-N (NPHASE = 2)
//   ViscousShockTube           viscous shock tube (Daru and Tenaud)     5Eq, 6Eq, 6Eq-N (NPHASE = 2)
//   AdvectionViscousDrop       viscous and conductive advection of a    6Eq, 6Eq-N3 (NPHASE = 3,
//                              drop                                     phase 3 absent)
//   AdvectionViscousMaterials  viscous and conductive advection of four 5Eq-N (NPHASE = 5)
//                              shapes in five materials
// =======================================
enum DiffusionProblem : int {
    Couette2Layer_ID             = 0,
    HeatConduction2Mat_ID        = 1,
    ViscousShockTube_ID          = 2,
    AdvectionViscousDrop_ID      = 3,
    AdvectionViscousMaterials_ID = 4
};

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS || PHYSICS == SIXEQS_IE_NPHASE)
#define DIFFUSION_HAS_COUETTE2LAYER
#endif
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS || (PHYSICS == SIXEQS_IE_NPHASE && NPHASE == 2))
#define DIFFUSION_HAS_HEATCONDUCTION2MAT
#define DIFFUSION_HAS_VISCOUSSHOCKTUBE
#endif
#if (PHYSICS == SIXEQS || PHYSICS == SIXEQS_IE_NPHASE)
#define DIFFUSION_HAS_ADVECTIONVISCOUSDROP
#endif
#if (PHYSICS == FIVEEQS_NPHASE && NPHASE == 5)
#define DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
#endif

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
#if (PHYSICS == SIXEQS)
    PhaseField_Parameter PhaseField_Parm{.ID_FluxPhase2Flux = 0};
#else
    PhaseField_Parameter PhaseField_Parm;
#endif

    //
    // Linear solver parameters (Grad_QPF)
    //
    LinearSystem_Parameter LinearSystem_Parm{.max_iter = 10'000, .max_fmg_iter = 10'000, .fixed_iter = 10'000, .verbose = 1};

    //
    // Equation of State (EOS): one value per phase, prob.EOS.gamma_1, prob.EOS.gamma_2, ... in
    // the two-phase models, the lists prob.EOS.gamma, ... in the N-phase models
    //
    int EOS = 1;

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
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
#else
    // Convex EOS fit parameters, and the EOS and transport parameters of each phase
    amrex::Array<amrex::Real, NPHASE> A, B, C, D;
    amrex::Array<amrex::Real, NPHASE> gamma, pinf, gGamma, pref, eref, cv, mu, muB, kappa;
#endif

    //
    // Governing equation parameters
    //
    amrex::Real K = 0.0;
    amrex::Array<amrex::Real, NSTATE> UpperBound;
    amrex::Array<amrex::Real, NSTATE> LowerBound;
    amrex::Array<int,        NSTATE> ID_BoundType;
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    amrex::Array<amrex::Real, 2>     p_min;
#else
    amrex::Array<amrex::Real, NPHASE>     p_min;
#endif

    //
    // Problem parameters (prob.* in the inputs file). Only the problem chosen by prob.problem
    // reads its parameters; a member shared by several problems (P0, eps, rho_1, rho_2) takes
    // that problem's default before its keys are read.
    //
    int problem = Couette2Layer_ID;

    // Shared: uniform pressure, trace volume fraction, phase densities
    amrex::Real P0 = 1.0;
    amrex::Real eps = 1.0e-8;
    amrex::Real rho_1 = 1.0, rho_2 = 10.0;

    // Couette2Layer: channel y_lo < y < y_hi, periodic in x, phase 1 below y = h and phase 2
    // above, the wall y = y_hi moving at U; prob.rho in the N-phase model. Set from the inputs:
    // the geometry and the exact steady shear stress tau = U/((h - y_lo)/mu_1 + (y_hi - h)/mu_2).
    amrex::Real U = 0.01;
    amrex::Real h = 0.5;
#if (PHYSICS == SIXEQS_IE_NPHASE || PHYSICS == FIVEEQS_NPHASE)
    amrex::Array<amrex::Real, NPHASE> rho_phase;
#endif
    amrex::Real y_lo = 0.0, y_hi = 1.0, L_x = 1.0, tau = 0.0;

    // HeatConduction2Mat: two materials in contact at x0 at rest, at T_L and T_R; the errors
    // are measured in the window |x - x0| < window. Set from the inputs: the cross-section area
    // and the exact solution (densities, cp, diffusivities chi, effusivities e, interface T_i).
    amrex::Real T_L = 1.005, T_R = 0.995;
    amrex::Real x0 = 0.0;
    amrex::Real window = 0.1875;
    amrex::Real Area = 1.0;
    amrex::Real rho_L = 1.0, rho_R = 1.0, cp_L = 1.0, cp_R = 1.0, chi_L = 1.0, chi_R = 1.0;
    amrex::Real e_L = 1.0, e_R = 1.0, T_i = 1.0;

    // AdvectionViscousDrop: a drop of phase 2, radius R, center (xc, yc), in phase 1, at the
    // uniform pressure P0 and velocity (u0, v0); densities rho_1 = 0.1 and rho_2 = 1 by default
    amrex::Real R = 0.4, xc = 0.0, yc = 0.0;
    amrex::Real u0 = 1.0, v0 = 1.0;

    // ViscousShockTube and AdvectionViscousMaterials have no prob.* parameters.

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

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
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
#else
   std::vector<amrex::Real> gammaV(NPHASE), pinfV(NPHASE), gGammaV(NPHASE), prefV(NPHASE), erefV(NPHASE), cvV(NPHASE), muV(NPHASE), muBV(NPHASE), kappaV(NPHASE);

   // The N-phase builds of this case differ in NPHASE (6Eq-N and 6Eq-N3): an inputs file run on
   // the wrong one gives a list of the wrong length
   auto queryphases = [&pp](const char* key, std::vector<amrex::Real>& values)
   {
      pp.queryarr(key, values);
      if (int(values.size()) != NPHASE){
         amrex::Abort(std::string("prob.EOS.") + key + " has " + std::to_string(values.size())
                      + " values, but this build has NPHASE = " + std::to_string(NPHASE)
                      + "; see the BUILD of this inputs file in variants.txt.");
      }
   };

   pp.query("EOS",EOS);
   queryphases("gamma",gammaV);
   queryphases("pinf",pinfV);

   queryphases("gGamma",gGammaV);
   queryphases("pref",prefV);
   queryphases("eref",erefV);

   queryphases("cv",cvV);

   queryphases("mu",muV);
   queryphases("muB",muBV);
   queryphases("kappa",kappaV);

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
#endif

   // ---------------------------------------
   // The problem
   // ---------------------------------------
   amrex::ParmParse pp_prob("prob");
   std::string name;
   pp_prob.get("problem",name);
   bool available = false;
   if (name == "Couette2Layer"){
      problem = Couette2Layer_ID;
#ifdef DIFFUSION_HAS_COUETTE2LAYER
      available = true;
#endif
   } else if (name == "HeatConduction2Mat"){
      problem = HeatConduction2Mat_ID;
#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
      available = true;
#endif
   } else if (name == "ViscousShockTube"){
      problem = ViscousShockTube_ID;
#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
      available = true;
#endif
   } else if (name == "AdvectionViscousDrop"){
      problem = AdvectionViscousDrop_ID;
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
      available = true;
#endif
   } else if (name == "AdvectionViscousMaterials"){
      problem = AdvectionViscousMaterials_ID;
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
      available = true;
#endif
   } else {
      amrex::Abort("prob.problem = " + name + " is not a problem of this case: Couette2Layer, "
                   "HeatConduction2Mat, ViscousShockTube, AdvectionViscousDrop or "
                   "AdvectionViscousMaterials.");
   }
   if (!available){
      amrex::Abort("prob.problem = " + name + " is not available in this build (model); see the "
                   "list of problems in prob/Parm.H and the BUILD of each inputs file in variants.txt.");
   }

   // ---------------------------------------
   // Its parameters
   // ---------------------------------------
#ifdef DIFFUSION_HAS_COUETTE2LAYER
   if (problem == Couette2Layer_ID)
   {
      pp_prob.query("U",U);
      pp_prob.query("h",h);
      pp_prob.query("P0",P0);
#if (PHYSICS == SIXEQS_IE_NPHASE)
      pp_prob.query("eps",eps);
      std::vector<amrex::Real> rhoV(NPHASE, 1.0);
      pp_prob.queryarr("rho",rhoV);
      if (NPHASE < 2 || int(rhoV.size()) != NPHASE){
         amrex::Abort("Couette2Layer needs NPHASE >= 2 and NPHASE values of prob.rho.");
      }
      for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
         rho_phase[ID_Phase] = rhoV[ID_Phase];
      }
#else
      pp_prob.query("rho_1",rho_1);
      pp_prob.query("rho_2",rho_2);
      pp_prob.query("eps",eps);
#endif

      amrex::ParmParse pp_geom("geometry");
      std::vector<amrex::Real> lo, hi;
      pp_geom.getarr("prob_lo",lo);
      pp_geom.getarr("prob_hi",hi);
      y_lo = lo[1];
      y_hi = hi[1];
      L_x  = hi[0] - lo[0];
      if (h <= y_lo || h >= y_hi){
          amrex::Abort("prob.h must lie inside the channel.");
      }
#if (PHYSICS == SIXEQS_IE_NPHASE)
      tau = U/((h - y_lo)/mu[0] + (y_hi - h)/mu[1]);
#else
      tau = U/((h - y_lo)/mu_1 + (y_hi - h)/mu_2);
#endif

      Initialize();
   }
#endif

#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
   if (problem == HeatConduction2Mat_ID)
   {
      pp_prob.query("P0",P0);
      pp_prob.query("T_L",T_L);
      pp_prob.query("T_R",T_R);
      pp_prob.query("x0",x0);
      pp_prob.query("eps",eps);
      pp_prob.query("window",window);

      amrex::ParmParse pp_geom("geometry");
      std::vector<amrex::Real> lo, hi;
      pp_geom.getarr("prob_lo",lo);
      pp_geom.getarr("prob_hi",hi);
      for (int d = 1; d < AMREX_SPACEDIM; ++d){ Area *= hi[d] - lo[d]; }

      Initialize();

      // Exact solution, after Initialize has set the EOS parameters A and B
#if (PHYSICS == SIXEQS_IE_NPHASE)
      rho_L = (A[0]*P0 + A[0]*B[0]/(1.0 + A[0]))/(cv[0]*T_L);
      rho_R = (A[1]*P0 + A[1]*B[1]/(1.0 + A[1]))/(cv[1]*T_R);
      cp_L  = cv[0]*(1.0 + A[0])/A[0];
      cp_R  = cv[1]*(1.0 + A[1])/A[1];
      chi_L = kappa[0]/(rho_L*cp_L);
      chi_R = kappa[1]/(rho_R*cp_R);
      e_L   = std::sqrt(kappa[0]*rho_L*cp_L);
      e_R   = std::sqrt(kappa[1]*rho_R*cp_R);
#else
      rho_L = (A1*P0 + A1*B1/(1.0 + A1))/(cv_1*T_L);
      rho_R = (A2*P0 + A2*B2/(1.0 + A2))/(cv_2*T_R);
      cp_L  = cv_1*(1.0 + A1)/A1;
      cp_R  = cv_2*(1.0 + A2)/A2;
      chi_L = kappa_1/(rho_L*cp_L);
      chi_R = kappa_2/(rho_R*cp_R);
      e_L   = std::sqrt(kappa_1*rho_L*cp_L);
      e_R   = std::sqrt(kappa_2*rho_R*cp_R);
#endif
      T_i   = (e_L*T_L + e_R*T_R)/(e_L + e_R);
   }
#endif

#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
   if (problem == ViscousShockTube_ID)
   {
      Initialize();
   }
#endif

#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
   if (problem == AdvectionViscousDrop_ID)
   {
      rho_1 = 0.1;
      rho_2 = 1.0;
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

#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
   if (problem == AdvectionViscousMaterials_ID)
   {
      Initialize();
   }
#endif
}

#endif
