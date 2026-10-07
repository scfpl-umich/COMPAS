// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef MACROS_H
#define MACROS_H

#include<AMReX.H>

using namespace amrex;

//==================================================================================================================================================================================================================
// Basic macro definitions
//==================================================================================================================================================================================================================
#define POW1(X) ((X))
#define POW2(X) ((X)*(X))
#define POW3(X) ((X)*(X)*(X))
#define TERNARY_SELECT_AMREX_SPACEDIM(X1, X2, X3) ((AMREX_SPACEDIM) == 1 ? (X1) : ((AMREX_SPACEDIM == 2) ? (X2) : (X3)))
#define POW_AMREX_SPACEDIM(X) TERNARY_SELECT_AMREX_SPACEDIM((POW1(X)),(POW2(X)),(POW3(X)))

// Create an alias for amrex::Real to clarify code
// #define scalar amrex::Real

/** Define a macro for the signum operator */
#define SGN(x) ((x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0))

#ifndef FIXED_DT
#define FIXED_DT false
#endif

#ifndef CONVERGENCE
#define CONVERGENCE false
#endif

#ifndef IMMERSED_BOUNDARY
#define IMMERSED_BOUNDARY false
#endif

#ifndef OUTPUT_TOOLS
#define OUTPUT_TOOLS false
#endif

#define PI 3.141592653589793238462643383279502884197
// #define MU0 1.25663706127e-6

#define FACES TERNARY_SELECT_AMREX_SPACEDIM(2,4,6)




//==================================================================================================================================================================================================================
// Physics macro definitions
//==================================================================================================================================================================================================================
#define FIVEEQS               (555)
#define SIXEQS                (666)
#define FIVEEQS_NPHASE        (55501)
#define SIXEQS_IE_NPHASE      (66601)


#ifndef PHYSICS
#define PHYSICS (FIVEEQS)
#endif

// Set some variables for the solver


// Five equations model, default diffusion to true
#if (PHYSICS == FIVEEQS)
// std::string PhysicsName = "FiveEquations";
#define NPHASE 1
#define NSTATE (int)(4 + AMREX_SPACEDIM)
#ifndef ADVECTION
#define ADVECTION (true)
#endif
#ifndef DIFFUSION
#define DIFFUSION (false)
#endif
#ifndef NONCONSERVATIVE
#define NONCONSERVATIVE (true)
#endif
#ifndef SURFACE_TENSION
#define SURFACE_TENSION (false)
#endif
// Non-conservative face component 0 holds u_n. With SURFACE_TENSION, the last AMREX_SPACEDIM
// components hold the capillary stress T.e_d (SurfaceTension_Parameter.H)
#if (SURFACE_TENSION == true)
#define NC_TERMS (1 + AMREX_SPACEDIM)
#else
#define NC_TERMS 1
#endif

#endif



// Five equations model with Phase-Fields and N-Phases
#if (PHYSICS == FIVEEQS_NPHASE)
// std::string PhysicsName = "FiveEquations";
#ifndef NPHASE
#define NPHASE 2
#endif
#define NSTATE (int)(2*NPHASE + AMREX_SPACEDIM + 1)
#ifndef ADVECTION
#define ADVECTION (true)
#endif
#ifndef DIFFUSION
#define DIFFUSION (false)
#endif
#ifndef NONCONSERVATIVE
#define NONCONSERVATIVE (true)
#endif
#ifndef SURFACE_TENSION
#define SURFACE_TENSION (false)
#endif
// Non-conservative face component 0 holds u_n. With SURFACE_TENSION, the last AMREX_SPACEDIM
// components hold the capillary stress T.e_d (SurfaceTension_Parameter.H)
#if (SURFACE_TENSION == true)
#define NC_TERMS (1 + AMREX_SPACEDIM)
#else
#define NC_TERMS 1
#endif

#endif


// Six equations model
#if (PHYSICS == SIXEQS)
// std::string PhysicsName = "SixEquations";
#define NPHASE 1
#define NSTATE (int)(5 + AMREX_SPACEDIM)
#ifndef ADVECTION
#define ADVECTION (true)
#endif
#ifndef DIFFUSION
#define DIFFUSION (false)
#endif
#ifndef NONCONSERVATIVE
#define NONCONSERVATIVE (true)
#endif
#ifndef SURFACE_TENSION
#define SURFACE_TENSION (false)
#endif
// Non-conservative face components 0 to 3: u_n, (u alpha1 p1)_n, (u alpha2 p2)_n and 0 in the
// hyperbolic step, and J1 rho1, J2 rho2, J1 rho1 EK, J2 rho2 EK in the Phase-Field step. With
// DIFFUSION, component 4 (INDEX_NC_ViscousWork) holds the viscous energy flux (u.tau)_n. With
// SURFACE_TENSION, the last AMREX_SPACEDIM + 1 components hold the capillary stress T.e_d and the
// capillary energy flux (u.T)_n (SurfaceTension_Parameter.H).
#if (DIFFUSION == true) && (SURFACE_TENSION == true)
#define NC_TERMS (6 + AMREX_SPACEDIM)
#elif (DIFFUSION == true)
#define NC_TERMS 5
#elif (SURFACE_TENSION == true)
#define NC_TERMS (5 + AMREX_SPACEDIM)
#else
#define NC_TERMS 4
#endif

#endif


// Six equations n-phase model (internal energy formulation)
#if (PHYSICS == SIXEQS_IE_NPHASE)
// std::string PhysicsName = "SixEquationsNphase";
#ifndef NPHASE
#define NPHASE 2
#endif
#define NSTATE (int)(3*NPHASE + AMREX_SPACEDIM + 1)
#ifndef ADVECTION
#define ADVECTION (true)
#endif
#ifndef DIFFUSION
#define DIFFUSION (false)
#endif
#ifndef NONCONSERVATIVE
#define NONCONSERVATIVE (true)
#endif
#ifndef SURFACE_TENSION
#define SURFACE_TENSION (false)
#endif
// Non-conservative face component 0 holds u_n. With SURFACE_TENSION, the last AMREX_SPACEDIM
// components hold the capillary stress T.e_d (SurfaceTension_Parameter.H)
#if (SURFACE_TENSION == true)
#define NC_TERMS (1 + AMREX_SPACEDIM)
#else
#define NC_TERMS 1
#endif

#endif

#endif

















