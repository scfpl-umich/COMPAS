// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SolveABecLaplacianAllLevels_H_
#define SolveABecLaplacianAllLevels_H_

inline
void SolveABecLaplacianAllLevels(amrex::Vector<amrex::MultiFab>& solution,
                                 amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& flux,
                                 amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& gradient,
                                 amrex::Vector<amrex::Geometry> const& geom,
                                 amrex::Vector<amrex::BoxArray> const& grids,
                                 amrex::Vector<amrex::DistributionMapping> const& dmap,
                                 amrex::Vector<amrex::IntVect> const& ref_ratio,
                                 amrex::Real const& ascalar,
                                 amrex::Real const& bscalar,
                                 amrex::Vector<amrex::MultiFab> const& acoef,
                                 amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> const& bcoef,
                                 amrex::Vector<amrex::MultiFab> const& rhs,
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
    bool composite_solve = LinearSystem_Parm.composite_solve;
    //Solve all levels together
    if (composite_solve)
    {
        SolveABecLaplacianComposite(solution, flux, gradient,
                                    geom, grids, dmap,
                                    ascalar, bscalar, acoef, bcoef, rhs,
                                    LinOpBC_lo, LinOpBC_hi,
                                    LevelBCData, RobinBC_a, RobinBC_b, RobinBC_f,
                                    ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                                    LinearSystem_Parm);
    }
    //Solve level by level
    else
    {
        const int nlevels = geom.size();
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            //Set coarse-fine BC
            if (ilev > 0) {
                ID_CoarseFineBC = 1;
                CoarseFineBCData = &(solution[ilev-1]);
                CoarseFineBCRatio = ref_ratio[ilev-1][0];
            }
            SolveABecLaplacianAtLevel(solution[ilev], flux[ilev], gradient[ilev],
                                      geom[ilev], grids[ilev], dmap[ilev],
                                      ascalar, bscalar, acoef[ilev], bcoef[ilev], rhs[ilev],
                                      LinOpBC_lo, LinOpBC_hi,
                                      LevelBCData[ilev], RobinBC_a[ilev], RobinBC_b[ilev], RobinBC_f[ilev],
                                      ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                                      LinearSystem_Parm);
        }
    }
    
}

#endif