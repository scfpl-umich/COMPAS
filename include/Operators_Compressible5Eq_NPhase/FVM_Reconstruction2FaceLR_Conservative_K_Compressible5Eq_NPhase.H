// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_K_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2FaceLR_Conservative_K_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_K(int i, int j, int k,
                                              amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeL,
                                              amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeR,
                                              amrex::Array4<amrex::Real const> const& Conservative,
                                              Parm const& EOS_Parameter,
                                              amrex::Real const& dx,
                                              int ID_Dimension,
                                              int ID_WENO,
                                              int ID_Limiter,
                                              int ID_Poly)
{
    //non-WENO
    if (ID_WENO == 0){
        //Godunov (1st order)
        if (ID_Limiter == 0){
            FVM_Reconstruction2FaceLR_Conservative_Godunov_K(i, j, k,
                                                             ConservativeL, ConservativeR,
                                                             Conservative,
                                                             EOS_Parameter,
                                                             ID_Dimension);
        }
        //MUSCL (2nd order)
        else{
            FVM_Reconstruction2FaceLR_Conservative_MUSCL_K(i, j, k,
                                                           ConservativeL, ConservativeR,
                                                           Conservative,
                                                           EOS_Parameter,
                                                           dx,
                                                           ID_Dimension, ID_Limiter);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            FVM_Reconstruction2FaceLR_Conservative_WENO3_K(i, j, k,
                                                           ConservativeL, ConservativeR,
                                                           Conservative,
                                                           EOS_Parameter,
                                                           dx,
                                                           ID_Dimension);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            FVM_Reconstruction2FaceLR_Conservative_WENO3B_K(i, j, k,
                                                            ConservativeL, ConservativeR,
                                                            Conservative,
                                                            EOS_Parameter,
                                                            dx,
                                                            ID_Dimension);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            FVM_Reconstruction2FaceLR_Conservative_WENO5_K(i, j, k,
                                                           ConservativeL, ConservativeR,
                                                           Conservative,
                                                           EOS_Parameter,
                                                           dx,
                                                           ID_Dimension);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            FVM_Reconstruction2FaceLR_Conservative_WENO5B_K(i, j, k,
                                                            ConservativeL, ConservativeR,
                                                            Conservative,
                                                            EOS_Parameter,
                                                            dx,
                                                            ID_Dimension, ID_Poly);
        }
    }
    
}

#endif
#endif