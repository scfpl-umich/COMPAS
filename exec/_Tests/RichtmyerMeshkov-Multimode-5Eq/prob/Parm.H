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

#include <Physics_Parameter_Compressible5Eq.H>
#include <FiniteVolume_Parameter.H>
#include <PhaseField_Parameter.H>
#include <LinearSystem_Parameter.H>

// Device storage for the interface modes. Their size is known only at run time, so they
// live outside Parm, which holds raw pointers to them.
namespace {
    amrex::Gpu::DeviceVector<amrex::Real> a_d;
    amrex::Gpu::DeviceVector<amrex::Real> kx_d;
    amrex::Gpu::DeviceVector<amrex::Real> phase_d;
}


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
    // Linear solver parameters
    //
    LinearSystem_Parameter LinearSystem_Parm{
        .max_iter = 10'000,
        .max_fmg_iter = 10'000,
        .fixed_iter = 10'000,
        .verbose = 1
    };

    //
    // Equation of State (EOS)
    //
    int EOS = 1;

    amrex::Real gamma_1, gamma_2;
    amrex::Real pinf_1, pinf_2;
    amrex::Real gGamma_1 = 0.0, gGamma_2 = 0.0, pref_1 = 0.0, pref_2 = 0.0, eref_1 = 0.0, eref_2 = 0.0;
    amrex::Real cv_1 = 1.0, cv_2 = 1.0;
    amrex::Real mu_1 = 0.0, mu_2 = 0.0, muB_1 = 0.0, muB_2 = 0.0, kappa_1 = 0.0, kappa_2 = 0.0;
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
    // Problem parameters (read from prob.* in the inputs file), nondimensional
    //
    amrex::Real rho1 = 1.0;    // light-gas (material 1) density ahead of the shock
    amrex::Real rho2 = 10.0;   // heavy-gas (material 2) density
    amrex::Real Mach = 1.5;    // shock Mach number in material 1
    amrex::Real a0   = 0.1;    // the shock starts at y = 1.1 a0

    // Interface modes, y = sum_m a[m] cos(kx[m] x + phase[m]), m < n_modes (device pointers)
    int n_modes = 0;
    const amrex::Real* a     = nullptr;
    const amrex::Real* kx    = nullptr;
    const amrex::Real* phase = nullptr;

    //
    // Initialization
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

    Initialize();

    amrex::ParmParse ppp("prob");

    ppp.query("rho1", rho1);
    ppp.query("rho2", rho2);
    ppp.query("Mach", Mach);
    ppp.query("a0",   a0);

    ppp.query("n_modes", n_modes);

    if (n_modes > 0)
    {
        amrex::Vector<amrex::Real> a_h, kx_h, phase_h;

        ppp.getarr("amplitudes",  a_h,  0, n_modes);
        ppp.getarr("wavenumbers", kx_h, 0, n_modes);

        if (ppp.contains("phases"))
            ppp.getarr("phases", phase_h, 0, n_modes);
        else
            phase_h.assign(n_modes, 0.0);

        a_d.resize(n_modes);
        kx_d.resize(n_modes);
        phase_d.resize(n_modes);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         a_h.begin(), a_h.end(),
                         a_d.begin());

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         kx_h.begin(), kx_h.end(),
                         kx_d.begin());

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         phase_h.begin(), phase_h.end(),
                         phase_d.begin());

        // Raw device pointers for the initial condition
        a     = a_d.data();
        kx    = kx_d.data();
        phase = phase_d.data();
    }
}


#endif
