// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SolveABecLaplacian_H_
#define SolveABecLaplacian_H_

inline
void SolveABecLaplacian(amrex::Vector<amrex::MultiFab*> const& solution,
                        amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& flux,
                        amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& gradient,
                        amrex::Vector<amrex::Geometry> const& geom,
                        amrex::Vector<amrex::BoxArray> const& grids,
                        amrex::Vector<amrex::DistributionMapping> const& dmap,
                        amrex::Real const& ascalar,
                        amrex::Real const& bscalar,
                        amrex::Vector<amrex::MultiFab const*> const& acoef,
                        amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM>> const& bcoef,
                        amrex::Vector<amrex::MultiFab const*> const& rhs,
                        amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_lo,
                        amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& LinOpBC_hi,
                        amrex::Vector<amrex::MultiFab const*> const& LevelBCData = {nullptr},
                        amrex::Vector<amrex::MultiFab const*> const& RobinBC_a = {nullptr},
                        amrex::Vector<amrex::MultiFab const*> const& RobinBC_b = {nullptr},
                        amrex::Vector<amrex::MultiFab const*> const& RobinBC_f = {nullptr},
                        int ID_CoarseFineBC = 0,
                        amrex::MultiFab const* CoarseFineBCData = nullptr,
                        int CoarseFineBCRatio = 2,
                        LinearSystem_Parameter const& LinearSystem_Parm = LinearSystem_Parameter{})
{
    //Construct MLABecLaplacian
    MLABecLaplacian mlabec;
    ConstructABecLaplacian(mlabec,
                           geom, grids, dmap,
                           ascalar, bscalar, acoef, bcoef,
                           LinOpBC_lo, LinOpBC_hi,
                           LevelBCData, RobinBC_a, RobinBC_b, RobinBC_f,
                           ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                           LinearSystem_Parm);

    //Linear solve
    Solver_MLMG(solution, flux, gradient,
                mlabec, rhs,
                LinearSystem_Parm);
    
}

#endif