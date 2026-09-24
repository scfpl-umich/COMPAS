// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxHyperbolic_K_Compressible5Eq_H_
#define FVM_Conservative2FluxHyperbolic_K_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxHyperbolic_K(int i, int j, int k,
                                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                                       amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
                                       amrex::Real& cmax,
                                       amrex::Array4<amrex::Real const> const& Conservative,
                                       Parm const& EOS_Parameter,
                                       amrex::Array<amrex::Real, AMREX_SPACEDIM> const& dx,
                                       int ID_Dimension,
                                       FiniteVolume_Parameter const& FiniteVolume_Parm,
                                       Physics_Parameter const& Physics_Parm)
{
    int ID_Bound   = FiniteVolume_Parm.ID_Bound;
    int ID_WENO    = FiniteVolume_Parm.ID_WENO;
    int ID_Limiter = FiniteVolume_Parm.ID_Limiter;
    int ID_Poly    = FiniteVolume_Parm.ID_Poly;
    int ID_THINC   = FiniteVolume_Parm.ID_THINC;
    int ID_Quad   = FiniteVolume_Parm.ID_Quad;

    int ID_K       = Physics_Parm.source_term;
    int ID_Riemann = Physics_Parm.ID_Riemann;


    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> ConservativeL;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> ConservativeR;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    if (ID_Bound == 0){
        if (ID_Quad){
            if (ID_THINC){
                amrex::Abort("[Error] THINC not implemented for WENO+Quadrature!\n");
            }
            if (!ID_WENO){
                amrex::Abort("[Error] Quadrature only defined for WENO5!\n");

            }
            FVM_Reconstruction2FaceLR_Conservative_WENO5_Quad_K(i, j, k,
                                                                Flux, Vel, cmax,
                                                                Conservative,
                                                                EOS_Parameter,
                                                                ID_Dimension,
                                                                ID_K, 
                                                                ID_Riemann);
        } else {
            if (ID_THINC && FiniteVolume_Parm.ID_THINC_Options){
                FVM_Reconstruction2Face_Conservative_THINCOptions_K(i, j, k,
                                                                    ConservativeL, ConservativeR,
                                                                    Conservative,
                                                                    EOS_Parameter,
                                                                    dx,
                                                                    ID_Dimension,
                                                                    FiniteVolume_Parm);
            } else if (ID_THINC){

                int ID_THINC_UseNormal   = FiniteVolume_Parm.ID_THINC_UseNormal;
                amrex::Real EPS = FiniteVolume_Parm.THINC_EPS;
                amrex::Real Beta = FiniteVolume_Parm.THINC_Beta;
                amrex::Real NormalThreshold = FiniteVolume_Parm.THINC_NormalThreshold;

                FVM_Reconstruction2Face_Conservative_THINC_K(i, j, k,
                                                            ConservativeL, ConservativeR,
                                                            Conservative,
                                                            EOS_Parameter,
                                                            dx,
                                                            ID_Dimension,
                                                            ID_THINC_UseNormal,
                                                            NormalThreshold, 
                                                            EPS, 
                                                            Beta,
                                                            ID_WENO, ID_Limiter, ID_Poly);
            } else {
                FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                        ConservativeL, ConservativeR,
                                                        Conservative,
                                                        EOS_Parameter,
                                                        dx[ID_Dimension],
                                                        ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);
            }

{BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::RiemannSolver");
            //Flux_Hyperbolic
            RiemannSolver(Flux, Vel, cmax,
                        ConservativeL, ConservativeR,
                        EOS_Parameter,
                        ID_Dimension, ID_K, ID_Riemann);
}
        }        
    }
    else if (ID_Bound == 1){
        if (ID_WENO == 1){
            FVM_Reconstruction2FaceLR_Conservative_WENO5BP_K(i, j, k,
                                                             ConservativeL, ConservativeR,
                                                             Conservative,
                                                             EOS_Parameter,
                                                             ID_Dimension, ID_Poly);
        }
        else{
            FVM_Reconstruction2FaceLR_Conservative_MUSCLBP_K(i, j, k,
                                                             ConservativeL, ConservativeR,
                                                             Conservative,
                                                             EOS_Parameter,
                                                             dx[ID_Dimension],
                                                             ID_Dimension, ID_Limiter);
        }

{BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::RiemannSolver");
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_m;
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_c;
        for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
            Conservative_c[i_State] = Conservative(i,j,k,i_State);
            //x
            if (ID_Dimension == 0){
                Conservative_m[i_State] = Conservative(i-1,j,k,i_State);
            }
            //y
            else if (ID_Dimension == 1){
                Conservative_m[i_State] = Conservative(i,j-1,k,i_State);
            }
            //z
            else if (ID_Dimension == 2){
                Conservative_m[i_State] = Conservative(i,j,k-1,i_State);
            }
        }

        //Flux_Hyperbolic
        int ID_Lambda = FiniteVolume_Parm.ID_Lambda;
        Riemann_HLL_Bound(Flux, Vel, cmax,
                          ConservativeL, ConservativeR,
                          Conservative_m, Conservative_c,
                          EOS_Parameter,
                          ID_Dimension, ID_K, ID_Lambda);

        amrex::Real Quadrature_Weight = FiniteVolume_Parm.Quadrature_Weight;
        amrex::Real CFLSpeed_c;
        Conservative2CFLSpeed_HLL(CFLSpeed_c,
                                  Conservative_c,
                                  cmax,
                                  Quadrature_Weight,
                                  EOS_Parameter,
                                  ID_K);
        amrex::Real CFLSpeed_m;
        Conservative2CFLSpeed_HLL(CFLSpeed_m,
                                  Conservative_m,
                                  cmax,
                                  Quadrature_Weight,
                                  EOS_Parameter,
                                  ID_K);
        cmax = std::max( CFLSpeed_c, CFLSpeed_m );
}
    }

    //FluxNC_Hyperbolic
    for (int iState = 0; iState < NC_TERMS; iState++){
        FluxNC[iState] = 0.0;
    }
    FluxNC[0] = Vel[ID_Dimension];
    
}







AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxHyperbolic_K(int i, int j, int k,
                                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                                       amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
                                       amrex::Real& cmax,
                                       amrex::Array4<amrex::Real const> const& Conservative,
                                       Parm const& EOS_Parameter,
                                       amrex::Real const& dx,
                                       int ID_Dimension,
                                       FiniteVolume_Parameter const& FiniteVolume_Parm,
                                       Physics_Parameter const& Physics_Parm)
{
    int ID_Bound   = FiniteVolume_Parm.ID_Bound;
    int ID_WENO    = FiniteVolume_Parm.ID_WENO;
    int ID_Limiter = FiniteVolume_Parm.ID_Limiter;
    int ID_Poly    = FiniteVolume_Parm.ID_Poly;

    int ID_K       = Physics_Parm.source_term;
    int ID_Riemann = Physics_Parm.ID_Riemann;


    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> ConservativeL;
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> ConservativeR;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    if (ID_Bound == 0){

        FVM_Reconstruction2FaceLR_Conservative_K(i, j, k,
                                                    ConservativeL, ConservativeR,
                                                    Conservative,
                                                    EOS_Parameter,
                                                    dx,
                                                    ID_Dimension, ID_WENO, ID_Limiter, ID_Poly);

{BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::RiemannSolver");
        //Flux_Hyperbolic
        RiemannSolver(Flux, Vel, cmax,
                      ConservativeL, ConservativeR,
                      EOS_Parameter,
                      ID_Dimension, ID_K, ID_Riemann);
}
    }
    else if (ID_Bound == 1){
        if (ID_WENO == 1){
            FVM_Reconstruction2FaceLR_Conservative_WENO5BP_K(i, j, k,
                                                             ConservativeL, ConservativeR,
                                                             Conservative,
                                                             EOS_Parameter,
                                                             ID_Dimension, ID_Poly);
        }
        else{
            FVM_Reconstruction2FaceLR_Conservative_MUSCLBP_K(i, j, k,
                                                             ConservativeL, ConservativeR,
                                                             Conservative,
                                                             EOS_Parameter,
                                                             dx,
                                                             ID_Dimension, ID_Limiter);
        }

{BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::RiemannSolver");
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_m;
        amrex::Array<amrex::Real,4+AMREX_SPACEDIM> Conservative_c;
        for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
            Conservative_c[i_State] = Conservative(i,j,k,i_State);
            //x
            if (ID_Dimension == 0){
                Conservative_m[i_State] = Conservative(i-1,j,k,i_State);
            }
            //y
            else if (ID_Dimension == 1){
                Conservative_m[i_State] = Conservative(i,j-1,k,i_State);
            }
            //z
            else if (ID_Dimension == 2){
                Conservative_m[i_State] = Conservative(i,j,k-1,i_State);
            }
        }

        //Flux_Hyperbolic
        int ID_Lambda = FiniteVolume_Parm.ID_Lambda;
        Riemann_HLL_Bound(Flux, Vel, cmax,
                          ConservativeL, ConservativeR,
                          Conservative_m, Conservative_c,
                          EOS_Parameter,
                          ID_Dimension, ID_K, ID_Lambda);

        amrex::Real Quadrature_Weight = FiniteVolume_Parm.Quadrature_Weight;
        amrex::Real CFLSpeed_c;
        Conservative2CFLSpeed_HLL(CFLSpeed_c,
                                  Conservative_c,
                                  cmax,
                                  Quadrature_Weight,
                                  EOS_Parameter,
                                  ID_K);
        amrex::Real CFLSpeed_m;
        Conservative2CFLSpeed_HLL(CFLSpeed_m,
                                  Conservative_m,
                                  cmax,
                                  Quadrature_Weight,
                                  EOS_Parameter,
                                  ID_K);
        cmax = std::max( CFLSpeed_c, CFLSpeed_m );
}
    }

    //FluxNC_Hyperbolic
    for (int iState = 0; iState < NC_TERMS; iState++){
        FluxNC[iState] = 0.0;
    }
    FluxNC[0] = Vel[ID_Dimension];
    
}

#endif
#endif