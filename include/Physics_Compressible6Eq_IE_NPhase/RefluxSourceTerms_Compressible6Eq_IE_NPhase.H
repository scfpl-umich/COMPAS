// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RefluxSourceTerms_Compressible6Eq_IE_NPhase_H_
#define RefluxSourceTerms_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void 
SaveSourceTermsBefore_Reflux (amrex::MultiFab & source_fab,
                             amrex::MultiFab const & dof, 
                             Parm const & EOS_Parameter)
{
    for (MFIter mfi(dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

    // ======== FLUX CALC AND UPDATE =========

        // Get the box and the box with a ghost layer
        const Box& bx = mfi.tilebox();

        // Pull the data into an array
        Array4<Real const> statein      = dof.const_array(mfi);
        Array4<Real      > source_arr   = source_fab.array(mfi);

        // Parm & EOS_Parameter = *d_parm;

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real c_avg = 1.0;
#if (SPACE == DGALERKIN)
            c_avg = 1./std::pow(std::sqrt(2.),amrex::Real(AMREX_SPACEDIM));
#endif
            amrex::Array<amrex::Real,NSTATE> Conservative;
            amrex::Array<amrex::Real,NSTATE> Primitive;

            for (int v = 0; v < NSTATE; v++){
                Conservative[v] = statein(i,j,k,v)*c_avg;
            }

            Conservative2Primitive(Primitive,Conservative,EOS_Parameter);

            amrex::Array<amrex::Real, NPHASE> Alpha;
            amrex::Array<amrex::Real, NPHASE> Alpha_rho;

            AMREX_D_TERM( amrex::Real u;, 
                        amrex::Real v;, 
                        amrex::Real w);

            amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
            amrex::Real P;

            Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), P, Primitive);


            for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
                source_arr(i,j,k,i_Phase) = -Alpha[i_Phase];
                source_arr(i,j,k,i_Phase+NPHASE) = Alpha_Ps[i_Phase];
            }
            
        });
    }


}



AMREX_GPU_HOST
AMREX_FORCE_INLINE
void 
RefluxSource_Terms (amrex::MultiFab & dof_new,
                   amrex::MultiFab & source_fab,
                   int lev,
                   amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr,
                   Geometry & geom,
                   BoxArray & ba,
                   DistributionMapping & dm,
                   Parm const & EOS_Parameter)
{

    /* Define the MF to store the velocity reflux */ 
    int ncomp = NC_TERMS; int ngrow = 0;
    MultiFab xi(ba,dm,ncomp,ngrow); xi.setVal(Real(0.0));

    /* Reflux the velocity flux */ 
    fr[lev+1]->Reflux(xi, 1.0, 0, 0, ncomp, geom);

    /* Reflux the phasic volume fractions and phasic pressures*/ 
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        amrex::MultiFab::Multiply(source_fab, xi, 0, i_Phase, 1, 0);    
        amrex::MultiFab::Saxpy(dof_new, Real(1.0), source_fab, i_Phase, i_Phase, 1, 0);

        amrex::MultiFab::Multiply(source_fab, xi, 0, i_Phase+NPHASE, 1, 0);    
        amrex::MultiFab::Saxpy(dof_new, Real(1.0), source_fab, i_Phase+NPHASE, i_Phase+2*NPHASE, 1, 0);
    }
   
}

#endif
#endif