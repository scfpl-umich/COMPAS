// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConstructABecLaplacian_H_
#define ConstructABecLaplacian_H_

inline
void ConstructABecLaplacian(amrex::MLABecLaplacian& mlabec,
                            amrex::Vector<amrex::Geometry> const& geom,
                            amrex::Vector<amrex::BoxArray> const& grids,
                            amrex::Vector<amrex::DistributionMapping> const& dmap,
                            amrex::Real const& ascalar,
                            amrex::Real const& bscalar,
                            amrex::Vector<amrex::MultiFab const*> const& acoef,
                            amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM>> const& bcoef,
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
    bool agglomeration           = LinearSystem_Parm.agglomeration;
    bool consolidation           = LinearSystem_Parm.consolidation;
    bool semicoarsening          = LinearSystem_Parm.semicoarsening;
    int max_coarsening_level     = LinearSystem_Parm.max_coarsening_level;
    int max_semicoarsening_level = LinearSystem_Parm.max_semicoarsening_level;

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

    const int nlevels = geom.size();

    mlabec.define(geom, grids, dmap, info);

    int LinOp_maxorder = LinearSystem_Parm.LinOp_maxorder;
    mlabec.setMaxOrder(LinOp_maxorder);

    //Set domain BC
    mlabec.setDomainBC(LinOpBC_lo,  //BC type at lower ends
                       LinOpBC_hi); //BC type at higher ends

    //Set coarse-fine BC
    if (ID_CoarseFineBC == 1) {
        mlabec.setCoarseFineBC(CoarseFineBCData, CoarseFineBCRatio);
    }

    //Set level BC
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        //for problem with inhomogeneous BC, we need to pass boundary values
        mlabec.setLevelBC(ilev, LevelBCData[ilev],
                          RobinBC_a[ilev], RobinBC_b[ilev], RobinBC_f[ilev]);
    }

    //Set scalars
    mlabec.setScalars(ascalar, bscalar);

    //Set coefficients
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        mlabec.setACoeffs(ilev, *(acoef[ilev]));
        mlabec.setBCoeffs(ilev, bcoef[ilev]);
    }
    
}

#endif