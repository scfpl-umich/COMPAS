// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SolveABecLaplacianComposite_HomBC_H_
#define SolveABecLaplacianComposite_HomBC_H_

inline
void SolveABecLaplacianComposite_HomBC(amrex::Vector<amrex::MultiFab>& solution,
                                       amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& flux,
                                       amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& gradient,
                                       amrex::Vector<amrex::Geometry> const& geom,
                                       amrex::Vector<amrex::BoxArray> const& grids,
                                       amrex::Vector<amrex::DistributionMapping> const& dmap,
                                       amrex::Real const& ascalar,
                                       amrex::Real const& bscalar,
                                       amrex::Vector<amrex::MultiFab> const& acoef,
                                       amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> const& bcoef,
                                       amrex::Vector<amrex::MultiFab> const& rhs,
                                       amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_lo,
                                       amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_hi,
                                       int ID_CoarseFineBC = 0,
                                       amrex::MultiFab const* CoarseFineBCData = nullptr,
                                       int CoarseFineBCRatio = 2,
                                       LinearSystem_Parameter const& LinearSystem_Parm = LinearSystem_Parameter{})
{
    //Construct MLABecLaplacian
    MLABecLaplacian mlabec;
    ConstructABecLaplacianAllLevels_HomBC(mlabec,
                                          geom, grids, dmap,
                                          ascalar, bscalar, acoef, bcoef,
                                          LinOpBC_lo, LinOpBC_hi,
                                          ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                                          LinearSystem_Parm);

    //Linear solve
    Solver_MLMG_AllLevels(solution, flux, gradient,
                          mlabec, rhs,
                          LinearSystem_Parm);
    
}

#endif