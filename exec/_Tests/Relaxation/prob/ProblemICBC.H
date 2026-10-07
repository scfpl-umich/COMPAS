// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

// The relaxation tests of the two builds of this case (BUILD in the GNUmakefile). The harnesses
// of the two six-equation models differ in most of their lines (the state layout, the
// relaxations and the checks), so each has its own file:
//   ProblemICBC_6Eq.H     two-phase six-equation model (SIXEQS)
//   ProblemICBC_6Eq-N.H   N-phase six-equation model (SIXEQS_IE_NPHASE, NPHASE = 3)
#include <Macros.H>
#if (PHYSICS == SIXEQS)
#include "ProblemICBC_6Eq.H"
#elif (PHYSICS == SIXEQS_IE_NPHASE)
#include "ProblemICBC_6Eq-N.H"
#endif

#endif  // IC_H_
