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

#include <Tools.H>
#include <STL_Reader.H>

// The STL queries are host only (see STL_Reader.H), and initial_condition calls
// them, so this case cannot be built for a GPU.
#ifdef AMREX_USE_GPU
#error "COMPAS-STL-5Eq is a CPU-only case: the STL queries in STL_Reader.H are host only."
#endif

inline StlGeom::STLObject stl;

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
    // Problem parameters (read from prob.* in the inputs file)
    //
    amrex::Real rho1 = 1.0;     // density of fluid 1, the light fluid above the pool
    amrex::Real rho2 = 100.0;   // density of fluid 2, the pool and the STL solid
    amrex::Real P0   = 100.0;   // pressure at y = 0
    AMREX_D_TERM(
    amrex::Real gx = 0.0;,      // gravity
    amrex::Real gy = -1.0;,
    amrex::Real gz = 0.0;)

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


    Initialize();

    // Problem parameters
    amrex::ParmParse ppp("prob");
    ppp.query("rho1", rho1);
    ppp.query("rho2", rho2);
    ppp.query("P0", P0);

    AMREX_D_TERM(
    ppp.query("gx", gx);,
    ppp.query("gy", gy);,
    ppp.query("gz", gz);)

    // --------------------------------------------------------------
    // STL solid: read the file, then place it. The mesh center is moved to
    // center_at, in the units of the file, and the mesh is then scaled about the
    // origin so that its largest dimension is normalize_size, and rotated by
    // rotate_about_point[3..5] degrees about x, y and z through rotate_about_point[0..2].
    // --------------------------------------------------------------
    amrex::ParmParse pps("prob.stl");
    std::string stl_filename;
    std::vector<amrex::Real> center_at(3, 0.0), rotate_about_point(6, 0.0);
    amrex::Real normalize_size = 1.0;
    pps.get("stl_filename", stl_filename);
    pps.queryarr("center_at", center_at, 0, 3);
    pps.queryarr("rotate_about_point", rotate_about_point, 0, 6);
    pps.query("normalize_size", normalize_size);

    stl.filename = stl_filename;
    if (!stl.ReadSTLFile()) {
        amrex::Abort("Could not read the STL file prob.stl.stl_filename = " + stl_filename);
    }

    stl.ResetTransform();
    stl.CenterAt({float(center_at[0]), float(center_at[1]), float(center_at[2])});
    stl.NormalizeSize(float(normalize_size));
    stl.RotateAroundPoint({float(rotate_about_point[0]), float(rotate_about_point[1]), float(rotate_about_point[2])},
                          {float(rotate_about_point[3]), float(rotate_about_point[4]), float(rotate_about_point[5])});

    // Bake the transforms into the vertices, then build the BVH the queries use.
    // stl.FitIntoBox({-0.4,-0.1,-0.1}, {0.4,0.1,0.1}) places the mesh in a given box
    // instead of the CenterAt / NormalizeSize pair above.
    stl.ApplyTransforms();
    stl.BuildBVH();

    // stl.PrintMeshSummary() reports the triangle count and the bounds, which is the
    // quickest way to see where your own geometry landed in the domain.
}

#endif
