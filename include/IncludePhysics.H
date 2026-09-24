// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef INCLUDE_PHYSICS_H_
#define INCLUDE_PHYSICS_H_

#include <Macros.H>

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
#include <Physics_Compressible_TwoPhase.H>
#include <Physics_Compressible_TwoPhase_PhaseField.H>

#if (PHYSICS == FIVEEQS)
#include <Physics_Compressible5Eq.H>
#include <Physics_Compressible5Eq_PhaseField.H>
#endif

#if (PHYSICS == SIXEQS)
#include <Physics_Compressible6Eq.H>
#include <Physics_Compressible6Eq_PhaseField.H>
#endif

#endif

#if (PHYSICS == FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)
#include <Physics_Compressible_NPhase.H>
#include <Physics_Compressible_NPhase_PhaseField.H>

#if (PHYSICS == FIVEEQS_NPHASE)
#include <Physics_Compressible5Eq_NPhase.H>
#include <Physics_Compressible5Eq_NPhase_PhaseField.H>
#endif

#if (PHYSICS == SIXEQS_IE_NPHASE)
#include <Physics_Compressible6Eq_IE_NPhase.H>
#include <Physics_Compressible6Eq_IE_NPhase_PhaseField.H>
#endif

#endif
#endif