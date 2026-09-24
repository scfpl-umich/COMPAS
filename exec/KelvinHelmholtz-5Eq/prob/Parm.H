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
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>

#include <Physics_Parameter_Compressible5Eq.H>
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
    // Initialization routine
    //
    void Initialize();
    #define DYNAMIC_INIT
    void DynamicInit();

    //
    // Problem parameters (read from prob.* in the inputs file), nondimensional
    //
    amrex::Real Mach   = 0.4;     // Blumen's Mach number, U/a (free-stream speed over the sound speed)
    amrex::Real U      = 1.0;     // free-stream speed, u = U tanh(y/delta)
    amrex::Real delta  = 1.0;     // half-thickness of the shear layer
    amrex::Real rho0   = 1.0;     // density of both materials
    amrex::Real eps    = 1.0e-5;  // amplitude of the seeded v, in units of U
    int         tracer = 1;       // 1: alpha_1 follows the layer; 0: material 1 everywhere
    amrex::Real tracer_width = 1.0; // tanh width of alpha_1 in units of delta; 1 follows the layer, dx/delta is sharp

    // Half-heights |y| < y_ref[lev] of the band refined from level lev to lev + 1
    amrex::Array<amrex::Real, 4> y_ref = {8.0, 6.0, 5.0, 4.0};

    // Derived in DynamicInit
    amrex::Real P0     = 1.0;     // uniform pressure, rho0 U^2 / (gamma M^2), so that a = U/M
    amrex::Real kx     = 1.0;     // wavenumber of the seeded mode, 2 pi / Lx (one wavelength in x)
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

    Initialize();

    amrex::ParmParse ppp("prob");
    ppp.query("Mach", Mach);
    ppp.query("U", U);
    ppp.query("delta", delta);
    ppp.query("rho0", rho0);
    ppp.query("eps", eps);
    ppp.query("tracer", tracer);
    ppp.query("tracer_width", tracer_width);
    amrex::Vector<amrex::Real> yr;
    ppp.queryarr("y_ref", yr);
    for (int l = 0; l < static_cast<int>(yr.size()) && l < 4; ++l) { y_ref[l] = yr[l]; }

    // The passive tracer needs the same fluid on both sides
    if (gamma_1 != gamma_2 || pinf_1 != pinf_2) {
        amrex::Abort("KelvinHelmholtz-5Eq: the two materials must have the same EOS, so that alpha_1 is a passive tracer.");
    }

    if (tracer_width <= 0.0) {
        amrex::Abort("KelvinHelmholtz-5Eq: prob.tracer_width must be positive.");
    }

    // Uniform pressure for the Mach number, and the wavenumber of one wavelength in x
    P0 = rho0 * U * U / (gamma_1 * Mach * Mach);
    amrex::ParmParse ppg("geometry");
    amrex::Vector<amrex::Real> lo(AMREX_SPACEDIM), hi(AMREX_SPACEDIM);
    ppg.getarr("prob_lo", lo, 0, AMREX_SPACEDIM);
    ppg.getarr("prob_hi", hi, 0, AMREX_SPACEDIM);
    kx = 2.0 * amrex::Math::pi<amrex::Real>() / (hi[0] - lo[0]);
    amrex::Print() << "Kelvin-Helmholtz: M = " << Mach << ", P0 = " << P0
                   << ", seeded wavenumber kx delta = " << kx * delta << "\n";
}

#endif
