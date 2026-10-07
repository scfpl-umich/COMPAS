// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VARIABLES_Compressible6Eq_H_
#define VARIABLES_Compressible6Eq_H_

// #include "Physics.H"
#include <Macros.H>
#include <Parm.H>
#include <Tools.H>

#if (PHYSICS == SIXEQS)

/* Define the variable names */
inline Vector<std::string> vnamesC = {"a1rho1","a2rho2",AMREX_D_DECL("rho_u", "rho_v", "rho_w"), "a1rho1E1", "a2rho2E2", "a1"};/**<Three-dimensional conservative variables*/
inline Vector<std::string> vnamesD = {"a1rho1","a2rho2","rho_u","rho_v","rho_w","a1rho1E1","a2rho2E2","a1","rho","u", "v", "w","P","a1P1","a2P2","T","c","Mach","vorticity_mag","q_criterion","lambda2"};/**<Available output variables*/

enum StateIndex{
    INDEX_Mass1           = 0,
    INDEX_Mass2           = 1,
    AMREX_D_DECL(
    INDEX_MomentumX       = 2,
    INDEX_MomentumY       = 3,
    INDEX_MomentumZ       = 4),
    INDEX_Energy1         = 2+AMREX_SPACEDIM,
    INDEX_Energy2         = 3+AMREX_SPACEDIM,
    INDEX_VolumeFraction1 = 4+AMREX_SPACEDIM
};/**< Enumerator for variable indices */

#if (DIFFUSION == true)
/** Non-conservative face component with the viscous energy flux (u.tau)_n, plus the heat flux of the
    mixture with the pressure-temperature relaxation (see NC_TERMS in Macros.H) */
enum NCIndex{
    INDEX_NC_ViscousWork  = 4
};
#endif

#endif
#endif