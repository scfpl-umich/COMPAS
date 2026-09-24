// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SolveABecLaplacianAtLevel_H_
#define SolveABecLaplacianAtLevel_H_

inline
void SolveABecLaplacianAtLevel(amrex::MultiFab& solution,
                               amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& flux,
                               amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& gradient,
                               amrex::Geometry const& geom,
                               amrex::BoxArray const& grids,
                               amrex::DistributionMapping const& dmap,
                               amrex::Real const& ascalar,
                               amrex::Real const& bscalar,
                               amrex::MultiFab const& acoef,
                               amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> const& bcoef,
                               amrex::MultiFab const& rhs,
                               amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_lo,
                               amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_hi,
                               amrex::MultiFab const* LevelBCData = nullptr,
                               amrex::MultiFab const* RobinBC_a = nullptr,
                               amrex::MultiFab const* RobinBC_b = nullptr,
                               amrex::MultiFab const* RobinBC_f = nullptr,
                               int ID_CoarseFineBC = 0,
                               amrex::MultiFab const* CoarseFineBCData = nullptr,
                               int CoarseFineBCRatio = 2,
                               LinearSystem_Parameter const& LinearSystem_Parm = LinearSystem_Parameter{})
{
    //Construct MLABecLaplacian
    MLABecLaplacian mlabec;
    ConstructABecLaplacianAtLevel(mlabec,
                                  geom, grids, dmap,
                                  ascalar, bscalar, acoef, bcoef,
                                  LinOpBC_lo, LinOpBC_hi,
                                  LevelBCData, RobinBC_a, RobinBC_b, RobinBC_f,
                                  ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                                  LinearSystem_Parm);

    //Linear solve
    Solver_MLMG_AtLevel(solution, flux, gradient,
                        mlabec, rhs,
                        LinearSystem_Parm);
    
}

#endif