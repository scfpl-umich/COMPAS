// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Conservative2FluxHyperbolic_K_Compressible6Eq_IE_NPhase_H_
#define FVM_Conservative2FluxHyperbolic_K_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Conservative2FluxHyperbolic_K(int i, int j, int k,
                                       amrex::Array<amrex::Real,NSTATE>& Flux,
                                       amrex::Array<amrex::Real,NC_TERMS>& FluxNC,
                                       amrex::Real& cmax,
                                       amrex::Array4<amrex::Real const> const& Conservative,
                                       Parm const& EOS_Parameter,
                                       amrex::Array<amrex::Real, AMREX_SPACEDIM> const& dx,
                                       int ID_Dimension,
                                       FiniteVolume_Parameter const& FiniteVolume_Parm,
                                       Physics_Parameter const& Physics_Parm)
{
    int ID_WENO    = FiniteVolume_Parm.ID_WENO;
    int ID_Limiter = FiniteVolume_Parm.ID_Limiter;
    int ID_Poly    = FiniteVolume_Parm.ID_Poly;
    int ID_THINC   = FiniteVolume_Parm.ID_THINC;

    amrex::Array<amrex::Real,NSTATE> ConservativeL;
    amrex::Array<amrex::Real,NSTATE> ConservativeR;

    if (ID_THINC){

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

    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    
    {
    BL_PROFILE("FVM_Conservative2FluxHyperbolic_D()::RiemannSolver");
    //Flux_Hyperbolic

    // for (int is=0; is<NSTATE; is++){
    //     amrex::Print()<<ConservativeL[is]<<" "<<ConservativeR[is]<<"\n";
    // }
    
    int ID_Riemann = Physics_Parm.ID_Riemann;
    RiemannSolver(Flux, Vel, cmax,
                  ConservativeL, ConservativeR,
                  EOS_Parameter,
                  ID_Dimension, ID_Riemann);
    }

    //FluxNC_Hyperbolic
    for (int iState = 0; iState < NC_TERMS; iState++){
        FluxNC[iState] = 0.0;
    }
    FluxNC[0] = Vel[ID_Dimension];
    
}

#endif
#endif