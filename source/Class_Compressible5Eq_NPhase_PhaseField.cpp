// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <COMPAS.H>
#include <Kernels.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Gpu.H>
#include <Tools.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <LinearAlgebra.H>
#include <Class_Compressible5Eq_NPhase_PhaseField.H>

#include <Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)
#include <IncludePhysics.H>
#include <Parm.H>
#include <Physics_Compressible_NPhase.H>
#include <Physics_Compressible5Eq_NPhase.H>
#include <Operators_Compressible5Eq_NPhase.H>
#include <Physics_PhaseField.H>
#include <Operators_PhaseField.H>
#include <Physics_Compressible_NPhase_PhaseField.H>
#include <Physics_Compressible5Eq_NPhase_PhaseField.H>
#include <Operators_Compressible5Eq_NPhase_PhaseField.H>
#endif

using namespace amrex;


Compressible5Eq_NPhase_PhaseField::~Compressible5Eq_NPhase_PhaseField ()
{
}



void 
Compressible5Eq_NPhase_PhaseField::SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
                                                                amrex::MultiFab const & dof)
{
#if (PHYSICS == FIVEEQS_NPHASE)
        Parm const* lparm = d_parm;

    for (MFIter mfi(dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // ======== FLUX CALC AND UPDATE =========

        // Get the box and the box with a ghost layer
        const Box& bx = mfi.tilebox();

        // Pull the data into an array
        Array4<Real const> statein      = dof.const_array(mfi);
        Array4<Real      > source_arr   = source_fab.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real c_avg = 1.0;
            amrex::Array<amrex::Real,NSTATE> Conservative;
            for (int iState = 0; iState < NSTATE; iState++){
                Conservative[iState] = statein(i,j,k,iState)*c_avg;
            }

            amrex::Array<amrex::Real,NPHASE> K;
            Conservative2K(K,
                           Conservative, *lparm,
                           lparm->Physics_Parm.source_term);

            for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
                amrex::Real Alphai = Conservative[INDEX_VolumeFraction1+i_Phase];
                source_arr(i,j,k,i_Phase) = Alphai*(1.0 + K[i_Phase]);
            }
            
        });
    }
#endif
}



void 
Compressible5Eq_NPhase_PhaseField::RefluxSourceTerms (amrex::MultiFab & new_dof,
                                                      amrex::MultiFab & source_fab,
                                                      int lev,
                                                      amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr)
{
#if (PHYSICS == FIVEEQS_NPHASE)
    int ncomp = 1; int ngrow = 0;

    const Geometry& geom_lev = geom[lev];
    const BoxArray ba_lev = grids[lev];
    const DistributionMapping dm_lev = dmap[lev];

    MultiFab xi(ba_lev,dm_lev,ncomp,ngrow); xi.setVal(Real(0.0));

    fr[lev+1]->Reflux(xi, 1.0, 0, 0, 1, geom_lev);

    /* Reflux the phasic volume fractions */ 
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::MultiFab::Multiply(source_fab, xi, 0, i_Phase, 1, 0);    
        amrex::MultiFab::Saxpy(new_dof, Real(-1.0), source_fab, i_Phase, INDEX_VolumeFraction1+i_Phase, 1, 0);
    }
#endif
}



void
Compressible5Eq_NPhase_PhaseField::RefluxLev (int lev)
{
#if (PHYSICS==FIVEEQS_NPHASE)
#if (NONCONSERVATIVE == true)
    MultiFab source_fab(grids[lev],dmap[lev],NPHASE,0);
    SaveSourceTermsBeforeReflux(source_fab, dof_new[lev]);
#endif
    // update lev based on coarse-fine flux mismatch
    flux_reg[lev+1]->Reflux(dof_new[lev], 1.0, 0, 0, NSTATE, geom[lev]);
            
#if (NONCONSERVATIVE == true)
    RefluxSourceTerms(dof_new[lev],source_fab,lev,flux_reg_nc);
#endif
#endif
}


void 
Compressible5Eq_NPhase_PhaseField::SetPhysicsBC ()
{
#if (PHYSICS==FIVEEQS_NPHASE)
    bcs.resize(NSTATE);

    //
    // Components are:
    //  Interior,        Inflow,          Outflow,          Symmetry,              SlipWall,            NoSlipWall
    //
    static int scalar_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
    };

    static int norm_vel_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd
    };

    static int tang_vel_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd
    };



    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            bcs[i_Phase].setLo(idim, scalar_bc[lo_bc[idim]]);    
            bcs[i_Phase].setHi(idim, scalar_bc[hi_bc[idim]]);

            bcs[NPHASE + AMREX_SPACEDIM + 1 + i_Phase].setLo(idim, scalar_bc[lo_bc[idim]]);
            bcs[NPHASE + AMREX_SPACEDIM + 1 + i_Phase].setHi(idim, scalar_bc[hi_bc[idim]]);

        }
    
        bcs[INDEX_Energy].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Energy].setHi(idim, scalar_bc[hi_bc[idim]]);
        
    }


#if (AMREX_SPACEDIM == 1)
    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);

#elif (AMREX_SPACEDIM == 2)
    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumX].setLo(1, tang_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumX].setHi(1, tang_vel_bc[hi_bc[1]]);

    bcs[INDEX_MomentumY].setLo(0, tang_vel_bc[lo_bc[0]]);    
    bcs[INDEX_MomentumY].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumY].setLo(1, norm_vel_bc[lo_bc[1]]);    
    bcs[INDEX_MomentumY].setHi(1, norm_vel_bc[hi_bc[1]]);

#elif (AMREX_SPACEDIM == 3)

    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumX].setLo(1, tang_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumX].setHi(1, tang_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumX].setLo(2, tang_vel_bc[lo_bc[2]]);
    bcs[INDEX_MomentumX].setHi(2, tang_vel_bc[hi_bc[2]]);

    bcs[INDEX_MomentumY].setLo(0, tang_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumY].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumY].setLo(1, norm_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumY].setHi(1, norm_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumY].setLo(2, tang_vel_bc[lo_bc[2]]);
    bcs[INDEX_MomentumY].setHi(2, tang_vel_bc[hi_bc[2]]);

    bcs[INDEX_MomentumZ].setLo(0, tang_vel_bc[lo_bc[0]]);    
    bcs[INDEX_MomentumZ].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumZ].setLo(1, tang_vel_bc[lo_bc[1]]);    
    bcs[INDEX_MomentumZ].setHi(1, tang_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumZ].setLo(2, norm_vel_bc[lo_bc[2]]);    
    bcs[INDEX_MomentumZ].setHi(2, norm_vel_bc[hi_bc[2]]);
#endif


#endif
}