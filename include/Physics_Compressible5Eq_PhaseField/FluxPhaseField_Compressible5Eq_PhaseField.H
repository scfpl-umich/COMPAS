// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Compressible5Eq_PhaseField_H_
#define FluxPhaseField_Compressible5Eq_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhaseField(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                    amrex::Real& cmax,
                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeL,
                    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeR,
                    amrex::Real const& Grad_Q,
                    Parm const& EOS_Parameter,
                    int ID_Flux,
                    int ID_FluxPhase2Flux)
{
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> ConservativePhaseL;
    Conservative2ConservativePhase(ConservativePhaseL,
                                   ConservativeL,
                                   EOS_Parameter);
    
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> ConservativePhaseR;
    Conservative2ConservativePhase(ConservativePhaseR,
                                   ConservativeR,
                                   EOS_Parameter);
    
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> FluxPhase;
    //UD
    if (ID_Flux == 0){
        FluxPhasePhaseField_UD(FluxPhase,
                               cmax,
                               ConservativePhaseL,
                               ConservativePhaseR,
                               Grad_Q,
                               EOS_Parameter);
    }
    //HLL
    else if (ID_Flux == 1){
        FluxPhasePhaseField_HLL(FluxPhase,
                                cmax,
                                ConservativePhaseL,
                                ConservativePhaseR,
                                Grad_Q,
                                EOS_Parameter);
    }
    
    //Phasic formulation
    if (ID_FluxPhase2Flux == 0){
        VectorPhase2Vector(Flux,
                           FluxPhase);
    }
    //Mixture formulation
    else if (ID_FluxPhase2Flux == 1){
        FluxPhasePhaseField2FluxPhaseField_Mixture(Flux,
                                                   FluxPhase,
                                                   ConservativeL,
                                                   ConservativeR,
                                                   EOS_Parameter);
    }
    
}

#endif