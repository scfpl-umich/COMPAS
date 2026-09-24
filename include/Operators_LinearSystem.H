// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Operators_LinearSystem_H_
#define Operators_LinearSystem_H_

#include<LinearSystem_Parameter.H>

#include<Solver_MLMG.H>
#include<Solver_MLMG_AtLevel.H>
#include<Solver_MLMG_AllLevels.H>

#include<AMReX_MLABecLaplacian.H>
#include<ConstructABecLaplacian.H>
#include<ConstructABecLaplacian_HomBC.H>

#include<ConstructABecLaplacianAtLevel.H>
#include<ConstructABecLaplacianAtLevel_HomBC.H>

#include<ConstructABecLaplacianAllLevels.H>
#include<ConstructABecLaplacianAllLevels_HomBC.H>

#include<SolveABecLaplacianAtLevel.H>
#include<SolveABecLaplacianAtLevel_HomBC.H>

#include<SolveABecLaplacianComposite.H>
#include<SolveABecLaplacianComposite_HomBC.H>

#include<SolveABecLaplacianAllLevels.H>
#include<SolveABecLaplacianAllLevels_HomBC.H>

#include<SolveABecLaplacian.H>
#include<SolveABecLaplacian_HomBC.H>

#include<SolveABecLaplacianLevels.H>
#include<SolveABecLaplacianLevels_HomBC.H>

#endif