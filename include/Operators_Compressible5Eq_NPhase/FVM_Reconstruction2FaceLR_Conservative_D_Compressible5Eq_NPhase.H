// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2FaceLR_Conservative_D_Compressible5Eq_NPhase_H_
#define FVM_Reconstruction2FaceLR_Conservative_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

template<typename T_Stencil>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2FaceLR_Conservative_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeL,
                                              amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& ConservativeR,
                                              T_Stencil const& Conservative_Stencil,
                                              Parm const& EOS_Parameter,
                                              amrex::Real const& dx,
                                              int ID_WENO,
                                              int ID_Limiter,
                                              int ID_Poly)
{
    //non-WENO
    if (ID_WENO == 0){
        //Godunov (1st order)
        if (ID_Limiter == 0){
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_m[i_State] = Conservative_Stencil[0][i_State];
                Conservative_c[i_State] = Conservative_Stencil[1][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_Godunov_D(ConservativeL, ConservativeR,
                                                             Conservative_m, Conservative_c,
                                                             EOS_Parameter);
        }
        //MUSCL (2nd order)
        else{
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_mm[i_State] = Conservative_Stencil[0][i_State];
                Conservative_m [i_State] = Conservative_Stencil[1][i_State];
                Conservative_c [i_State] = Conservative_Stencil[2][i_State];
                Conservative_p [i_State] = Conservative_Stencil[3][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_MUSCL_D(ConservativeL, ConservativeR,
                                                           Conservative_mm, Conservative_m, Conservative_c, Conservative_p,
                                                           EOS_Parameter,
                                                           dx,
                                                           ID_Limiter);
        }
    }
    //WENO
    else if (ID_WENO == 1){
        //WENO3 (3rd order)
        if (ID_Limiter == 3){
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_mm[i_State] = Conservative_Stencil[0][i_State];
                Conservative_m [i_State] = Conservative_Stencil[1][i_State];
                Conservative_c [i_State] = Conservative_Stencil[2][i_State];
                Conservative_p [i_State] = Conservative_Stencil[3][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_WENO3_D(ConservativeL, ConservativeR,
                                                           Conservative_mm, Conservative_m, Conservative_c, Conservative_p,
                                                           EOS_Parameter,
                                                           dx);
        }
        //WENO3B (3rd order Bound)
        else if (ID_Limiter == 4){
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_mm[i_State] = Conservative_Stencil[0][i_State];
                Conservative_m [i_State] = Conservative_Stencil[1][i_State];
                Conservative_c [i_State] = Conservative_Stencil[2][i_State];
                Conservative_p [i_State] = Conservative_Stencil[3][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_WENO3B_D(ConservativeL, ConservativeR,
                                                            Conservative_mm, Conservative_m, Conservative_c, Conservative_p,
                                                            EOS_Parameter,
                                                            dx);
        }
        //WENO5 (5th order)
        else if (ID_Limiter == 5){
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mmm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_pp;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_mmm[i_State] = Conservative_Stencil[0][i_State];
                Conservative_mm [i_State] = Conservative_Stencil[1][i_State];
                Conservative_m  [i_State] = Conservative_Stencil[2][i_State];
                Conservative_c  [i_State] = Conservative_Stencil[3][i_State];
                Conservative_p  [i_State] = Conservative_Stencil[4][i_State];
                Conservative_pp [i_State] = Conservative_Stencil[5][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_WENO5_D(ConservativeL, ConservativeR,
                                                           Conservative_mmm, Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp,
                                                           EOS_Parameter,
                                                           dx);
        }
        //WENO5B (5th order Bound)
        else if (ID_Limiter == 6){
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mmm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_mm;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_m;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_c;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_p;
            amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative_pp;
            for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
                Conservative_mmm[i_State] = Conservative_Stencil[0][i_State];
                Conservative_mm [i_State] = Conservative_Stencil[1][i_State];
                Conservative_m  [i_State] = Conservative_Stencil[2][i_State];
                Conservative_c  [i_State] = Conservative_Stencil[3][i_State];
                Conservative_p  [i_State] = Conservative_Stencil[4][i_State];
                Conservative_pp [i_State] = Conservative_Stencil[5][i_State];
            }
            FVM_Reconstruction2FaceLR_Conservative_WENO5B_D(ConservativeL, ConservativeR,
                                                            Conservative_mmm, Conservative_mm, Conservative_m, Conservative_c, Conservative_p, Conservative_pp,
                                                            EOS_Parameter,
                                                            dx,
                                                            ID_Poly);
        }
    }
    
}

#endif
#endif