// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseField_RefluxSourceTerms_Compressible6Eq_PhaseField_H_
#define PhaseField_RefluxSourceTerms_Compressible6Eq_PhaseField_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)


// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// PhaseField_SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
//                                         amrex::MultiFab const & dof, 
//                                         Parm * EOS_Parameter,
//                                         PHYS * p)
// {
//     MultiFab::Copy(source_fab, dof, 0, 0, NSTATE, 0);
// }




// template <typename PHYS>
// AMREX_GPU_HOST
// AMREX_FORCE_INLINE
// void 
// PhaseField_RefluxSourceTerms (amrex::MultiFab & dof_new,
//                               amrex::MultiFab const & source_fab,
//                               int lev,
//                               amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr,
//                               Geometry & geom,
//                               BoxArray & ba,
//                               DistributionMapping & dm,
//                               Parm & EOS_Parameter,
//                               PHYS * p)
// {
//     int ncomp = NC_TERMS; int ngrow = 0;

//     MultiFab xi(ba,dm,ncomp,ngrow); xi.setVal(Real(0.0));

//     fr[lev+1]->Reflux(xi, 1.0, 0, 0, ncomp, geom);

// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//     {
//         for (MFIter mfi(dof_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();
            
//             // Pull the data into an array
//             Array4<Real const> ConservativeArray = source_fab.const_array(mfi);
//             Array4<Real const> xiArray           = xi.const_array(mfi);
//             Array4<Real      > dofArray          = dof_new.array(mfi);

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
//                                        EOS_Parameter);
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
// }

#endif
#endif