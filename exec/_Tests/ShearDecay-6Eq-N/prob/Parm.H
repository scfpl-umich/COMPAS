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
    // Shear wave u = U sin(k (x_a - a_lo)), k = 2 pi mode/L_a, with x_a the coordinate along
    // the axis prob.axis (1: y, or 2: z in 3D), in a uniform mixture of the phases at rest
    // otherwise, at the uniform pressure P0
    amrex::Real U = 0.01;
    int mode = 1;
    int axis = 1;                        // Direction of the wave vector: 1 (y) or, in 3D, 2 (z)
    amrex::Real P0 = 1.0;
    amrex::Array<amrex::Real, NPHASE> alpha_phase; // Volume fractions (prob.alpha), normalized to sum to 1
    amrex::Array<amrex::Real, NPHASE> rho_phase;   // Phase densities (prob.rho)

    // Set from geometry.prob_lo and geometry.prob_hi
    amrex::Real a_lo = 0.0, L_a = 1.0, Volume = 1.0;

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
   pp_prob.query("U",U);
   pp_prob.query("mode",mode);
   pp_prob.query("axis",axis);
   if (axis < 1 || axis >= AMREX_SPACEDIM){
       amrex::Abort("prob.axis = '" + std::to_string(axis) + "' is not recognized. Valid values are 1 (y) and, in 3D, 2 (z).");
   }
   pp_prob.query("P0",P0);
   std::vector<amrex::Real> alphaV(NPHASE, 1.0/NPHASE), rhoV(NPHASE, 1.0);
   pp_prob.queryarr("alpha",alphaV);
   pp_prob.queryarr("rho",rhoV);
   if (int(alphaV.size()) != NPHASE || int(rhoV.size()) != NPHASE){
      amrex::Abort("prob.alpha and prob.rho need NPHASE values each.");
   }
   amrex::Real Sum_Alpha = 0.0;
   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){ Sum_Alpha += alphaV[ID_Phase]; }
   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
      alpha_phase[ID_Phase] = alphaV[ID_Phase]/Sum_Alpha;
      rho_phase[ID_Phase]   = rhoV[ID_Phase];
   }

   amrex::ParmParse pp_geom("geometry");
   std::vector<amrex::Real> lo, hi;
   pp_geom.getarr("prob_lo",lo);
   pp_geom.getarr("prob_hi",hi);
   a_lo = lo[axis];
   L_a  = hi[axis] - lo[axis];
   Volume = 1.0;
   for (int d = 0; d < AMREX_SPACEDIM; ++d){
       Volume *= hi[d] - lo[d];
   }

   Initialize();

}

#endif
