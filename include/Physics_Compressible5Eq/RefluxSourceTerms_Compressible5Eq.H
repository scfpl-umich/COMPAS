// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RefluxSourceTerms_Compressible5Eq_H_
#define RefluxSourceTerms_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)


// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
//                              amrex::MultiFab const & dof, 
//                              Parm * EOS_Parameter,
//                              PHYS * p)
// {

//     int const source = p->source_term;

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
//             amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative;
//             for (int iState = 0; iState < NSTATE; iState++){
//                 Conservative[iState] = statein(i,j,k,iState);
//             }
//             amrex::Real Alpha1 = Conservative[INDEX_VolumeFraction1];
//             amrex::Array<amrex::Real,2> K;
//             Conservative2K(K,
//                            Conservative, *EOS_Parameter,
//                            source);
//             amrex::Real K1 = K[0];
//             source_arr(i,j,k,0) = Alpha1*(1.0 + K1);
//         });
//     }
// }




// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// RefluxSourceTerms (amrex::MultiFab & dof_new,
//                    amrex::MultiFab const & source_fab,
//                    int lev,
//                    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr,
//                    Geometry & geom,
//                    BoxArray & ba,
//                    DistributionMapping & dm,
//                    Parm & EOS_Parameter,
//                    PHYS * p)
// {
//     int ncomp = 1; int ngrow = 0;

//     MultiFab xi(ba,dm,ncomp,ngrow); xi.setVal(Real(0.0));

//     fr[lev+1]->Reflux(xi, 1.0, 0, 0, 1, geom);

//     // Update source term here via xi 
//     amrex::MultiFab::Multiply(xi, source_fab, 0, 0, 1, 0);

//     amrex::MultiFab::Saxpy(dof_new, Real(-1.0), xi, 0, INDEX_VolumeFraction1, 1, 0);
// }









#endif
#endif