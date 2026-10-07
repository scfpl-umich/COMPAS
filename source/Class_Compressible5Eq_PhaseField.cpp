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
#include <LinearAlgebra.H>
#include <Tools.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <Class_Compressible5Eq_PhaseField.H>

#include <Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==FIVEEQS_NPHASE)
#include <IncludePhysics.H>
#include <Parm.H>
#include <Physics_PhaseField.H>
#include <Operators_PhaseField.H>

#if (PHYSICS==FIVEEQS)
#include <Physics_Compressible_TwoPhase.H>
#include <Physics_Compressible_TwoPhase_PhaseField.H>
#endif

#if (PHYSICS==FIVEEQS_NPHASE)
#include <Physics_Compressible_NPhase.H>
#include <Physics_Compressible_NPhase_PhaseField.H>
#endif

#endif

#if (PHYSICS==FIVEEQS)
#include <Physics_Compressible5Eq.H>
#include <Operators_Compressible5Eq.H>
#include <Physics_Compressible5Eq_PhaseField.H>
#include <Operators_Compressible5Eq_PhaseField.H>
#endif

#if (PHYSICS==FIVEEQS_NPHASE)
#include <Physics_Compressible5Eq_NPhase.H>
#include <Operators_Compressible5Eq_NPhase.H>
#include <Physics_Compressible5Eq_NPhase_PhaseField.H>
#include <Operators_Compressible5Eq_NPhase_PhaseField.H>
#endif


using namespace amrex;


Compressible5Eq_PhaseField::~Compressible5Eq_PhaseField ()
{
}

void
Compressible5Eq_PhaseField::TimeStepChecksDerived(int step)
{
#if (PHYSICS == FIVEEQS || PHYSICS == FIVEEQS_NPHASE) 
    if (h_parm->Physics_Parm.source_term != 1 &&
        step >= h_parm->Physics_Parm.source_term_step_switch && t_new[0] >= h_parm->Physics_Parm.source_term_time_switch){
        h_parm->Physics_Parm.source_term = 1;
        // The kernels read the device copy of the parameters
        Parm* lparm = d_parm;
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
        {
            lparm->Physics_Parm.source_term = 1;
        });
    }
#endif
}
    

void 
Compressible5Eq_PhaseField::SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
                             amrex::MultiFab const & dof)
{
#if (PHYSICS == FIVEEQS)
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
            amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative;
            for (int iState = 0; iState < NSTATE; iState++){
                Conservative[iState] = statein(i,j,k,iState);
            }
            amrex::Real Alpha1 = Conservative[INDEX_VolumeFraction1];
            // Kapila coefficient alpha_1 (1 + K_1) = phi_1 (1 + K_1) for NASG (see NonConservative)
            amrex::Real Phi1 = (lparm->Physics_Parm.source_term == 1)
                             ? NASG_Phi(Alpha1, Conservative[INDEX_Mass1], lparm->C1) : Alpha1;
            amrex::Array<amrex::Real,2> K;
            Conservative2K(K,
                           Conservative, *lparm,
                           lparm->Physics_Parm.source_term);
            amrex::Real K1 = K[0];
            source_arr(i,j,k,0) = Phi1*(1.0 + K1);
        });
    }
#endif
}



void 
Compressible5Eq_PhaseField::RefluxSourceTerms (amrex::MultiFab & new_dof,
                                               amrex::MultiFab & source_fab,
                                               int lev,
                                               amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr)
{
#if (PHYSICS == FIVEEQS)
    int ncomp = 1; int ngrow = 0;

    const Geometry& geom_lev = geom[lev];
    const BoxArray ba_lev = grids[lev];
    const DistributionMapping dm_lev = dmap[lev];

    MultiFab xi(ba_lev,dm_lev,ncomp,ngrow); xi.setVal(Real(0.0));

    fr[lev+1]->Reflux(xi, 1.0, 0, 0, 1, geom_lev);

    // Update source term here via xi 
    amrex::MultiFab::Multiply(xi, source_fab, 0, 0, 1, 0);

    amrex::MultiFab::Saxpy(new_dof, Real(-1.0), xi, 0, INDEX_VolumeFraction1, 1, 0);
#endif
}


void
Compressible5Eq_PhaseField::RefluxLev (int lev)
{
#if (PHYSICS==FIVEEQS)
#if (NONCONSERVATIVE == true)
    MultiFab source_fab(grids[lev],dmap[lev],NC_TERMS,0);
    SaveSourceTermsBeforeReflux(source_fab, dof_new[lev]);
#endif
#if (SURFACE_TENSION == true)
    MultiFab state_before(grids[lev],dmap[lev],NSTATE,0);
    MultiFab::Copy(state_before, dof_new[lev], 0, 0, NSTATE, 0);
#endif
    // update lev based on coarse-fine flux mismatch
    flux_reg[lev+1]->Reflux(dof_new[lev], 1.0, 0, 0, NSTATE, geom[lev]);
            
#if (NONCONSERVATIVE == true)
    RefluxSourceTerms(dof_new[lev],source_fab,lev,flux_reg_nc);
#endif
#if (SURFACE_TENSION == true)
    SurfaceTension_RefluxSourceTerms(dof_new[lev], state_before, lev, flux_reg_nc);
#endif
#endif
}


void 
Compressible5Eq_PhaseField::SetPhysicsBC ()
{
#if (PHYSICS==FIVEEQS)
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
        bcs[INDEX_Mass1].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Mass2].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Energy].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_VolumeFraction1].setLo(idim, scalar_bc[lo_bc[idim]]);

        bcs[INDEX_Mass1].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_Mass2].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_Energy].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_VolumeFraction1].setHi(idim, scalar_bc[hi_bc[idim]]);
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