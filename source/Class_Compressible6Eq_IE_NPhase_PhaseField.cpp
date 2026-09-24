// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.



#include <COMPAS.H>
#include <Kernels.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Gpu.H>
#include "Tools.H"
#include "Operators_FV.H"
#include "Operators_FD.H"
#include <LinearAlgebra.H>
#include <Class_Compressible6Eq_IE_NPhase_PhaseField.H>

#include <Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)
#include "Parm.H"
#include <IncludePhysics.H>
// #include "Physics_Compressible_TwoPhase.H"
#include "Physics_Compressible6Eq_IE_NPhase.H"
#include <Operators_Compressible6Eq_IE_NPhase.H>

// #include <Physics_Compressible_TwoPhase_PhaseField.H>
#include <Physics_Compressible6Eq_IE_NPhase_PhaseField.H>
#include <Operators_Compressible6Eq_IE_NPhase_PhaseField.H>
#endif

using namespace amrex;


Compressible6Eq_IE_NPhase_PhaseField::~Compressible6Eq_IE_NPhase_PhaseField ()
{
}

void
Compressible6Eq_IE_NPhase_PhaseField::TimeStepChecksDerived(int step)
{  
}

void 
Compressible6Eq_IE_NPhase_PhaseField::PostTimeStage (MultiFab & mf, 
                                           int step, 
                                           amrex::Real time, 
                                           int lev)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
    // Correction_Energy(mf, d_parm);
    // The flag is read on the host from h_parm. Relaxation_Pressure copies *parm by value
    // on the host into its kernel, so it takes h_parm, not d_parm
    if (h_parm->Physics_Parm.pressure_relaxation){
        Relaxation_Pressure(mf, h_parm);
    }
#endif
}

void 
Compressible6Eq_IE_NPhase_PhaseField::PhaseField_PostTimeStage (MultiFab & mf, 
                                                      int step, 
                                                      amrex::Real time, 
                                                      int lev)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
    PostTimeStage(mf, step, time, lev);
#endif
}

void 
Compressible6Eq_IE_NPhase_PhaseField::SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
                                                         amrex::MultiFab const & dof)
{
// #if (PHYSICS == SIXEQS_IE_NPHASE)
//     // MultiFab::Copy(source_fab, dof, 0, 0, NSTATE, 0);
//     Parm const* lparm = d_parm;

//     SaveSourceTermsBefore_Reflux (source_fab,
//                                  dof, 
//                                  *lparm);
// #endif

#if (PHYSICS == SIXEQS_IE_NPHASE)
    MultiFab::Copy(source_fab, dof, 0, 0, NSTATE, 0);
#endif

}



void 
Compressible6Eq_IE_NPhase_PhaseField::RefluxSourceTerms (amrex::MultiFab & new_dof,
                                               amrex::MultiFab & source_fab,
                                               int lev,
                                               amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr)
{

// #if (PHYSICS == SIXEQS_IE_NPHASE)
// Parm const* lparm = d_parm;

// RefluxSource_Terms ( new_dof,
//                     source_fab,
//                     lev,
//                     fr,
//                     geom[lev],
//                     grids[lev],
//                     dmap[lev],
//                     *lparm);
// #endif

#if (PHYSICS == SIXEQS_IE_NPHASE)
    int ncomp = NC_TERMS; int ngrow = 0;

    Parm const* EOS_Parameter = d_parm;

    const Geometry& geom_lev = geom[lev];
    const BoxArray ba_lev = grids[lev];
    const DistributionMapping dm_lev = dmap[lev];

    MultiFab xi(ba_lev,dm_lev,ncomp,ngrow); xi.setVal(Real(0.0));

    fr[lev+1]->Reflux(xi, 1.0, 0, 0, ncomp, geom_lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(new_dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            
            // Pull the data into an array
            Array4<Real const> ConservativeArray = source_fab.const_array(mfi);
            Array4<Real const> xiArray           = xi.const_array(mfi);
            Array4<Real      > dofArray          = new_dof.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Array<Real,NSTATE> Conservative;
                for (int iState = 0; iState < NSTATE; ++iState){
                    Conservative[iState] = ConservativeArray(i,j,k,iState);
                }
                Array<Real,NSTATE> Primitive;
                Conservative2Primitive(Primitive,
                                       Conservative,
                                       *EOS_Parameter);

                amrex::Array<amrex::Real,NPHASE> Alpha;
                amrex::Array<amrex::Real,NPHASE> Alpha_rho;
                amrex::Array<amrex::Real,NPHASE> Alpha_Ps;
                AMREX_D_TERM( amrex::Real u;,
                              amrex::Real v;,
                              amrex::Real w );
                amrex::Real p;

                amrex::Array<amrex::Real,NPHASE> Alpha_rho_e;

                AMREX_D_TERM( amrex::Real rho_u;, 
                              amrex::Real rho_v;, 
                              amrex::Real rho_w);
                            
                amrex::Real rho_E;

                Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u, v, w), p, Primitive);

                amrex::Real Div_Vel  = xiArray(i,j,k,0);

                for (int i_Phase=0; i_Phase<NPHASE; ++i_Phase){
                    dofArray(i,j,k,a1+i_Phase) += -Alpha[i_Phase]*Div_Vel;
                    dofArray(i,j,k,a1r1e1+i_Phase) += Alpha_Ps[i_Phase]*Div_Vel;
                }
            });
        } // end mfi
    } // end omp
#endif
}

void
Compressible6Eq_IE_NPhase_PhaseField::RefluxLev (int lev)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
#if (NONCONSERVATIVE == true)
    MultiFab source_fab(grids[lev],dmap[lev],NSTATE,0);
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
Compressible6Eq_IE_NPhase_PhaseField::PhaseField_SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
                                                                    amrex::MultiFab const & dof)
{
#if (PHYSICS == SIXEQS_IE_NPHASE)
    // MultiFab::Copy(source_fab, dof, 0, 0, NSTATE, 0);
#endif
}



void 
Compressible6Eq_IE_NPhase_PhaseField::PhaseField_RefluxSourceTerms (amrex::MultiFab & new_dof,
                                                          amrex::MultiFab const & source_fab,
                                                          int lev,
                                                          amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr)
{
// #if (PHYSICS == SIXEQS_IE_NPHASE)
//     int ncomp = NC_TERMS; int ngrow = 0;

//     Parm const* EOS_Parameter = d_parm;

//     const Geometry& geom_lev = geom[lev];
//     const BoxArray ba_lev = grids[lev];
//     const DistributionMapping dm_lev = dmap[lev];

//     MultiFab xi(ba_lev,dm_lev,ncomp,ngrow); xi.setVal(Real(0.0));

//     fr[lev+1]->Reflux(xi, 1.0, 0, 0, ncomp, geom_lev);

// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         for (MFIter mfi(new_dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();
            
//             // Pull the data into an array
//             Array4<Real const> ConservativeArray = source_fab.const_array(mfi);
//             Array4<Real const> xiArray           = xi.const_array(mfi);
//             Array4<Real      > dofArray          = new_dof.array(mfi);

//             amrex::ParallelFor(bx,
//             [=] AMREX_GPU_DEVICE (int i, int j, int k)
//             {
//                 Array<Real,NSTATE> Conservative;
//                 for (int iState = 0; iState < NSTATE; ++iState){
//                     Conservative[iState] = ConservativeArray(i,j,k,iState);
//                 }
//                 Array<Real,NSTATE> Primitive;
//                 Conservative2Primitive(Primitive,
//                                        Conservative,
//                                        *EOS_Parameter);
//                 amrex::Real Alpha1_rho1;
//                 amrex::Real Alpha2_rho2;
//    AMREX_D_TERM(amrex::Real u;,
//                 amrex::Real v;,
//                 amrex::Real w);
//                 amrex::Real Alpha1_p1;
//                 amrex::Real Alpha2_p2;
//                 amrex::Real Alpha1;
//                 Primitive2State(Alpha1_rho1, Alpha2_rho2,
//                    AMREX_D_DECL(u, v, w),
//                                 Alpha1_p1, Alpha2_p2,
//                                 Alpha1,
//                                 Primitive);

//                 amrex::Real rho = Alpha1_rho1 + Alpha2_rho2;
//                 amrex::Real EK = 0.5*(AMREX_D_TERM(u*u, + v*v, + w*w));

//                 amrex::Real Div_J1_rho1     = -xiArray(i,j,k,0);
//                 amrex::Real Div_J2_rho2     = -xiArray(i,j,k,1);
//                 amrex::Real Div_J1_rho1_EK1 = -xiArray(i,j,k,2);
//                 amrex::Real Div_J2_rho2_EK2 = -xiArray(i,j,k,3);

//                 amrex::Real Vel_dot_J1_rho1_dot_Grad_Vel = Div_J1_rho1_EK1 - EK*Div_J1_rho1;
//                 amrex::Real Vel_dot_J2_rho2_dot_Grad_Vel = Div_J2_rho2_EK2 - EK*Div_J2_rho2;

//                 dofArray(i,j,k,INDEX_Energy1) +=  (Alpha1_rho1/rho*Vel_dot_J2_rho2_dot_Grad_Vel - Alpha2_rho2/rho*Vel_dot_J1_rho1_dot_Grad_Vel);
//                 dofArray(i,j,k,INDEX_Energy2) += -(Alpha1_rho1/rho*Vel_dot_J2_rho2_dot_Grad_Vel - Alpha2_rho2/rho*Vel_dot_J1_rho1_dot_Grad_Vel);
//             });
//         } // end mfi
//     } // end omp
// #endif
}

void
Compressible6Eq_IE_NPhase_PhaseField::PhaseField_RefluxLev (int lev)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
#if (NONCONSERVATIVE == true)
    MultiFab source_fab(grids[lev],dmap[lev],NSTATE,0);
    PhaseField_SaveSourceTermsBeforeReflux(source_fab, dof_new[lev]);
#endif
    // update lev based on coarse-fine flux mismatch
    flux_reg[lev+1]->Reflux(dof_new[lev], 1.0, 0, 0, NSTATE, geom[lev]);

#if (NONCONSERVATIVE == true)
    PhaseField_RefluxSourceTerms(dof_new[lev],source_fab,lev,flux_reg_nc);
#endif
#endif
}


void 
Compressible6Eq_IE_NPhase_PhaseField::SetPhysicsBC ()
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
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


    for(int iState=0; iState < 3*NPHASE; ++iState){
        for (int idim=0; idim < AMREX_SPACEDIM; ++idim){
            bcs[iState].setLo(idim, scalar_bc[lo_bc[idim]]);
            bcs[iState].setHi(idim, scalar_bc[hi_bc[idim]]);
        }
    }

    for (int idim=0; idim < AMREX_SPACEDIM; ++idim){
        bcs[INDEX_Total_Energy].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Total_Energy].setHi(idim, scalar_bc[hi_bc[idim]]);
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


void 
Compressible6Eq_IE_NPhase_PhaseField::Correction_Energy (MultiFab& mf,
			                                   Parm const* parm)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
			// Get the box
            const Box& bx = mfi.tilebox();

            // Pull the data into an array
            Array4<Real const> statein   = mf.const_array(mfi);
            Array4<Real      > stateout  = mf.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
				Array<Real,NSTATE> Conservative;
				for (int iState = 0; iState < NSTATE; ++iState){
					Conservative[iState] = statein(i,j,k,iState);
				}

				Array<Real,NSTATE> ConservativeC;
				EnergyCorrection(ConservativeC,
					             Conservative,
								 *parm);

				for (int iState = 0; iState < NSTATE; ++iState){
					stateout(i,j,k,iState) = ConservativeC[iState];
				}
            });
        } // end mfi
    } // end omp
#endif
}


void 
Compressible6Eq_IE_NPhase_PhaseField::Relaxation_Pressure (MultiFab& mf,
				                                 Parm const* parm)
{
#if (PHYSICS==SIXEQS_IE_NPHASE)

    PressureRelaxationInstantaneous ( mf, *parm);

#endif
}