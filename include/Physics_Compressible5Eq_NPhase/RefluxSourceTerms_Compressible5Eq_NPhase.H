// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RefluxSourceTerms_Compressible5Eq_NPhase_H_
#define RefluxSourceTerms_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
//                              amrex::MultiFab const & dof, 
//                              Parm & EOS_Parameter,
//                              PHYS * p)
// {
//     int source = p->source_term;

//     for (MFIter mfi(dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//     {
//         // ======== FLUX CALC AND UPDATE =========

//         // Get the box and the box with a ghost layer
//         const Box& bx = mfi.tilebox();

//         // Pull the data into an array
//         Array4<Real const> statein      = dof.const_array(mfi);
//         Array4<Real      > source_arr   = source_fab.array(mfi);

//         amrex::ParallelFor(bx,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k)
//         {
//             amrex::Real c_avg = 1.0;
//             amrex::Array<amrex::Real,NSTATE> Conservative;
//             for (int iState = 0; iState < NSTATE; iState++){
//                 Conservative[iState] = statein(i,j,k,iState)*c_avg;
//             }

//             amrex::Array<amrex::Real,NPHASE> K;
//             Conservative2K(K,
//                            Conservative, EOS_Parameter,
//                            source);

//             for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
//                 amrex::Real Alphai = Conservative[INDEX_VolumeFraction1+i_Phase];
//                 source_arr(i,j,k,i_Phase) = Alphai*(1.0 + K[i_Phase]);
//             }
            
//         });
//     }
// }




// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// RefluxSourceTerms (amrex::MultiFab & dof_new,
//                    amrex::MultiFab & source_fab,
//                    int lev,
//                    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr,
//                    Geometry & geom,
//                    BoxArray & ba,
//                    DistributionMapping & dm,
//                    Parm & EOS_Parameter,
//                    PHYS * p)
// {
//     /* Define the MF to store the velocity reflux */ 
//     int ncomp = 1; int ngrow = 0;
//     MultiFab xi(ba,dm,ncomp,ngrow); xi.setVal(Real(0.0));

//     /* Reflux the velocity flux */ 
//     fr[lev+1]->Reflux(xi, 1.0, 0, 0, 1, geom);

//     /* Reflux the phasic volume fractions */ 
//     for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
//         amrex::MultiFab::Multiply(source_fab, xi, 0, i_Phase, 1, 0);    
//         amrex::MultiFab::Saxpy(dof_new, Real(-1.0), source_fab, i_Phase, INDEX_VolumeFraction1+i_Phase, 1, 0);
//     }
    

// }









#endif
#endif