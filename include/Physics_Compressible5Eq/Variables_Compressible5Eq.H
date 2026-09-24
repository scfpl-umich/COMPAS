// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VARIABLES_Compressible5Eq_H_
#define VARIABLES_Compressible5Eq_H_

// #include "Physics.H"
#include <Macros.H>
#include <Parm.H>
#include <Tools.H>

#if (PHYSICS == FIVEEQS)

/* Define the variable names */
inline Vector<std::string> vnamesC = {"a1rho1","a2rho2",AMREX_D_DECL("rho_u", "rho_v", "rho_w"), "rho_E", "a1"};/**<Three-dimensional conservative variables*/
inline Vector<std::string> vnamesD = {"a1rho1","a2rho2","rho_u","rho_v","rho_w","rho_E","a1","rho", "u", "v", "w", "P","T","c","Mach","vorticity_mag","q_criterion","lambda2"};/**<Available output variables*/

enum StateIndex{
    INDEX_Mass1           = 0,
    INDEX_Mass2           = 1,
    AMREX_D_DECL(
    INDEX_MomentumX       = 2,
    INDEX_MomentumY       = 3,
    INDEX_MomentumZ       = 4),
    INDEX_Energy          = 2+AMREX_SPACEDIM,
    INDEX_VolumeFraction1 = 3+AMREX_SPACEDIM
};/**< Enumerator for variable indices */

#endif
#endif