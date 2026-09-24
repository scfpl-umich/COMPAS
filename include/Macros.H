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
#define NC_TERMS 1

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
#define NC_TERMS 1

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
#if (DIFFUSION == true)
#error "DIFFUSION=true is not implemented for PHYSICS=SIXEQS yet."
#endif
#define NC_TERMS 4

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
#if (DIFFUSION == true)
#error "DIFFUSION=true is not yet wired into the core for PHYSICS=SIXEQS_IE_NPHASE: the diffusion kernels in Operators_Compressible6Eq_IE_NPhase exist, but the viscous-work coupling (tau) into the non-conservative terms has not been ported."
#endif
#define NC_TERMS 1

#endif

#endif

















