// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>

// =======================================
// Viscous and conductive fluxes (DIFFUSION = true): five problems in one case, chosen at run
// time by prob.problem (Parm.H). Each problem is in its own header, in a namespace of its
// name, with the code of each model under #if (PHYSICS == ...). The functions below are the
// ones COMPAS calls; they call those of the problem of the run.
// =======================================

#include "Couette2Layer.H"
#include "HeatConduction2Mat.H"
#include "ViscousShockTube.H"
#include "AdvectionViscousDrop.H"
#include "AdvectionViscousMaterials.H"


// =======================================
// Initial condition
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
{
    switch (parm.problem)
    {
#ifdef DIFFUSION_HAS_COUETTE2LAYER
    case Couette2Layer_ID:
        Couette2Layer::initial_condition(Conservative, sharp, x, y, z,
                                         AMREX_D_DECL(dx, dy, dz),
                                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                         parm);
        break;
#endif
#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
    case HeatConduction2Mat_ID:
        HeatConduction2Mat::initial_condition(Conservative, sharp, x, y, z,
                                              AMREX_D_DECL(dx, dy, dz),
                                              AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                              parm);
        break;
#endif
#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
    case ViscousShockTube_ID:
        ViscousShockTube::initial_condition(Conservative, sharp, x, y, z,
                                            AMREX_D_DECL(dx, dy, dz),
                                            AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                            parm);
        break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
    case AdvectionViscousDrop_ID:
        AdvectionViscousDrop::initial_condition(Conservative, sharp, x, y, z,
                                                AMREX_D_DECL(dx, dy, dz),
                                                AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                                parm);
        break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
    case AdvectionViscousMaterials_ID:
        AdvectionViscousMaterials::initial_condition(Conservative, sharp, x, y, z,
                                                     AMREX_D_DECL(dx, dy, dz),
                                                     AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                                     parm);
        break;
#endif
    default:
        amrex::ignore_unused(Conservative, sharp, x, y, z, AMREX_D_DECL(dx, dy, dz),
                             AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));
        amrex::Abort("Diffusion: the problem is not available in this build");
    }
}


// =======================================
// Boundary condition: only Couette2Layer has a user condition (its moving no-slip walls)
// =======================================

#ifdef DIFFUSION_HAS_COUETTE2LAYER
#define USER_BOUNDARY_FUNC_ROBIN
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_Robin(int i, int j, int k, int n,
                            Real& BC_Dirichlet,
                            Real& BC_Neumann,
                            Real& BC_Function,
                            Real const& xBC,
                            Real const& yBC,
                            Real const& zBC,
                            Real const& time,
                            Array4<Real const> const& q,
               AMREX_D_DECL(int ilo, int jlo, int klo),
               AMREX_D_DECL(int ihi, int jhi, int khi),
               AMREX_D_DECL(Real const& dx,
                            Real const& dy,
                            Real const& dz),
                            Parm const* parm)
{
    if (parm->problem == Couette2Layer_ID) {
        Couette2Layer::BoundaryFunction_Robin(i, j, k, n, BC_Dirichlet, BC_Neumann, BC_Function,
                                              xBC, yBC, zBC, time, q,
                                              AMREX_D_DECL(ilo, jlo, klo),
                                              AMREX_D_DECL(ihi, jhi, khi),
                                              AMREX_D_DECL(dx, dy, dz),
                                              parm);
        return;
    }
    amrex::ignore_unused(i, j, k, n, BC_Dirichlet, BC_Neumann, BC_Function, xBC, yBC, zBC, time, q,
                         AMREX_D_DECL(ilo, jlo, klo), AMREX_D_DECL(ihi, jhi, khi), AMREX_D_DECL(dx, dy, dz));
    amrex::Abort("Diffusion: only Couette2Layer has a user boundary condition (run.lo_bc, run.hi_bc)");
}
#endif


// =======================================
// In-situ diagnostics: the CSV file of the problem (couette.csv, conduction.csv, box.csv,
// advection.csv or materials.csv)
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    switch (UserOutput::HostParm(o).problem)
    {
#ifdef DIFFUSION_HAS_COUETTE2LAYER
    case Couette2Layer_ID:             Couette2Layer::UserOutputFunction(o);             break;
#endif
#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
    case HeatConduction2Mat_ID:        HeatConduction2Mat::UserOutputFunction(o);        break;
#endif
#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
    case ViscousShockTube_ID:          ViscousShockTube::UserOutputFunction(o);          break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
    case AdvectionViscousDrop_ID:      AdvectionViscousDrop::UserOutputFunction(o);      break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
    case AdvectionViscousMaterials_ID: AdvectionViscousMaterials::UserOutputFunction(o); break;
#endif
    default: break;
    }
}

#endif  // IC_H_
