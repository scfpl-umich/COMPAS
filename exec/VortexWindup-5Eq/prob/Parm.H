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

#include <string>

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
    amrex::Real beta   = 5.0;     // vortex strength
    amrex::Real rc     = 1.0;     // core radius
    amrex::Real xc     = 0.0;     // center of the vortex
    amrex::Real yc     = 0.0;
    amrex::Real rho0   = 1.0;     // density far from the vortex (both materials)
    amrex::Real p0     = 1.0;     // pressure far from the vortex
    amrex::Real angle  = 0.0;     // direction of the normal n of the interface, in degrees from x
    amrex::Real offset = 0.0;     // distance of the interface from the center, along n
    amrex::Real width  = 0.0;     // 0: material 2 fills the half plane n.(x - c) > offset;
                                  // > 0: the slab 0 < n.(x - c) - offset < width
    amrex::Real alpha_min = 1.0e-8; // smallest volume fraction of either material
    int         nsub   = 4;       // samples per direction for the in-situ L1 error

    // Derived in DynamicInit
    amrex::Real nx = 1.0, ny = 0.0;   // n = (cos angle, sin angle)
    amrex::Real Omega0 = 0.0;         // angular velocity at the center, beta e^(1/2) / (2 pi rc)
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
    ppp.query("beta", beta);
    ppp.query("rc", rc);
    ppp.query("xc", xc);
    ppp.query("yc", yc);
    ppp.query("rho0", rho0);
    ppp.query("p0", p0);
    ppp.query("angle", angle);
    ppp.query("offset", offset);
    ppp.query("width", width);
    ppp.query("alpha_min", alpha_min);
    ppp.query("nsub", nsub);

    // The vortex is steady only if the two materials are the same ideal gas
    if (EOS != 0 || gamma_1 != gamma_2 || pinf_1 != 0.0 || pinf_2 != 0.0) {
        amrex::Abort("VortexWindup-5Eq: the two materials must be the same ideal gas (prob.EOS.EOS = 0, "
                     "gamma_1 = gamma_2, pinf = 0), so that the interface is passive.");
    }
    if (rc <= 0.0 || width < 0.0 || nsub < 1) {
        amrex::Abort("VortexWindup-5Eq: prob.rc must be positive, prob.width not negative and prob.nsub at least 1.");
    }
    // The core temperature must stay positive
    const amrex::Real Tmin = 1.0 - (gamma_1 - 1.0) * beta * beta * rho0 / (8.0 * gamma_1 * M_PI * M_PI * p0) * std::exp(1.0);
    if (Tmin <= 0.0) {
        amrex::Abort("VortexWindup-5Eq: prob.beta is too large for prob.p0/prob.rho0; the core temperature is negative.");
    }

    const amrex::Real ph = angle * M_PI / 180.0;
    nx = std::cos(ph);
    ny = std::sin(ph);
    Omega0 = beta * std::exp(0.5) / (2.0 * M_PI * rc);
    amrex::Print() << "Vortex wind-up: beta = " << beta << ", rc = " << rc << ", Omega(0) = " << Omega0
                   << " (one turn of the center in " << 2.0 * M_PI / Omega0 << "), T(0)/T_inf = " << Tmin << "\n";
}

#endif
