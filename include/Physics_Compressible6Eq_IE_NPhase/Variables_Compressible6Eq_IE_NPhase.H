// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VARIABLES_Compressible6Eq_IE_NPhase_H_
#define VARIABLES_Compressible6Eq_IE_NPhase_H_

// #include "Physics.H"
#include <Macros.H>
#include <Parm.H>
#include <Tools.H>

#if (PHYSICS == SIXEQS_IE_NPHASE)

/* Define the variable names */
inline Vector<std::string> vnamesC; /**<Three-dimensional conservative variables*/
inline Vector<std::string> vnamesD; /**<Available output variables*/

enum StateIndex{
    INDEX_VolumeFraction1   = 0,
    INDEX_Mass1             = NPHASE,
    INDEX_Energy1           = 2*NPHASE,
    AMREX_D_DECL(
    INDEX_MomentumX         = 3*NPHASE,
    INDEX_MomentumY         = 3*NPHASE+1,
    INDEX_MomentumZ         = 3*NPHASE+2),
    INDEX_Total_Energy      = 3*NPHASE+AMREX_SPACEDIM,
};/**< Enumerator for variable indices */

enum stateIndex {a1 = 0, a1r1 = NPHASE, a1r1e1 = 2*NPHASE,
                AMREX_D_DECL(ru=3*NPHASE, rv = 3*NPHASE+1, rw = 3*NPHASE+2), 
                rE = 3*NPHASE+ AMREX_SPACEDIM};

#endif
#endif