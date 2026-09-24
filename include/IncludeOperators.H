// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef _INCLUDE_OPERATORS_H_
#define _INCLUDE_OPERATORS_H_


#if (PHYSICS == FIVEEQS)
#include<Operators_Compressible5Eq.H>
#include<Operators_Compressible5Eq_PhaseField.H>
#endif

#if (PHYSICS == FIVEEQS_NPHASE)
#include<Operators_Compressible5Eq_NPhase.H>
#include<Operators_Compressible5Eq_NPhase_PhaseField.H>
#endif

#if (PHYSICS == SIXEQS)
#include<Operators_Compressible6Eq.H>
#include<Operators_Compressible6Eq_PhaseField.H>
#endif

#if (PHYSICS == SIXEQS_IE_NPHASE)
#include<Operators_Compressible6Eq_IE_NPhase.H>
#include<Operators_Compressible6Eq_IE_NPhase_PhaseField.H>
#endif




#endif