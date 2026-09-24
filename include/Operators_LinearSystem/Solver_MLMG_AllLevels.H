// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Solver_MLMG_AllLevels_H_
#define Solver_MLMG_AllLevels_H_

inline
void Solver_MLMG_AllLevels(amrex::Vector<amrex::MultiFab>& solution,
                           amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& flux,
                           amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& gradient,
                           amrex::MLLinOp& mllinop,
                           amrex::Vector<amrex::MultiFab> const& rhs,
                           LinearSystem_Parameter const& LinearSystem_Parm)
{
    int max_iter                      = LinearSystem_Parm.max_iter;
    int max_fmg_iter                  = LinearSystem_Parm.max_fmg_iter;
    int fixed_iter                    = LinearSystem_Parm.fixed_iter;
    int verbose                       = LinearSystem_Parm.verbose;
    int bottom_verbose                = LinearSystem_Parm.bottom_verbose;
    amrex::BottomSolver bottom_solver = LinearSystem_Parm.bottom_solver;

    MLMG mlmg(mllinop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setFixedIter(fixed_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setBottomSolver(bottom_solver);

    amrex::Real tol_rel = LinearSystem_Parm.tol_rel;
    amrex::Real tol_abs = LinearSystem_Parm.tol_abs;
    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    mlmg.getFluxes(GetVecOfArrOfPtrs(flux));
    mlmg.getGradSolution (GetVecOfArrOfPtrs(gradient));
    
}

#endif