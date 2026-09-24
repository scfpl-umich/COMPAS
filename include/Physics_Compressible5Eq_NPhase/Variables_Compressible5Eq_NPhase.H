// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VARIABLES_Compressible5Eq_NPhase_H_
#define VARIABLES_Compressible5Eq_NPhase_H_

// #include "Physics.H"
#include <Macros.H>
#include <LoopMacros.H>
#include <Parm.H>
#include <Tools.H>

#if (PHYSICS==FIVEEQS_NPHASE)

inline Vector<std::string> vnamesC;
inline Vector<std::string> vnamesD;

enum StateIndex{
    INDEX_Mass1           = 0,
    AMREX_D_DECL(
    INDEX_MomentumX       = NPHASE,
    INDEX_MomentumY       = NPHASE+1,
    INDEX_MomentumZ       = NPHASE+2),
    INDEX_Energy          = NPHASE+AMREX_SPACEDIM,
    INDEX_VolumeFraction1 = NPHASE+AMREX_SPACEDIM+1
};/**< Enumerator for variable indices */

#endif
#endif