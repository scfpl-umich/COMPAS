// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Compressible6Eq_IE_NPhase_PhaseField_H_
#define FluxPhaseField_Compressible6Eq_IE_NPhase_PhaseField_H_

#if (PHYSICS == SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhaseField(amrex::Array<amrex::Real, NSTATE>& Flux,
                    amrex::Real& cmax,
                    amrex::Array<amrex::Real, NSTATE> const& ConservativeL,
                    amrex::Array<amrex::Real, NSTATE> const& ConservativeR,
                    amrex::Array<amrex::Real,NPHASE> const& Grad_Q,
                    Parm const& EOS_Parameter,
                    int ID_FluxPhaseField,
                    int ID_FluxPhase2Flux,
                    int ID_BE)
{
     amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> ConservativePhaseL;
    Conservative2ConservativePhase(ConservativePhaseL,
                                   ConservativeL,
                                   EOS_Parameter);
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> ConservativePhaseR;
    Conservative2ConservativePhase(ConservativePhaseR,
                                   ConservativeR,
                                   EOS_Parameter);
    
    amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> FluxPhase;
    //UD
    if (ID_FluxPhaseField == 0){
        //Base formulation
        if (ID_BE == 1){
            FluxPhasePhaseField_UD_Base(FluxPhase,
                                        cmax,
                                        ConservativePhaseL,
                                        ConservativePhaseR,
                                        Grad_Q,
                                        EOS_Parameter);
        }
        //Efficient formulation
        else if (ID_BE == 2){
            FluxPhasePhaseField_UD_Efficient(FluxPhase,
                                             cmax,
                                             ConservativePhaseL,
                                             ConservativePhaseR,
                                             Grad_Q,
                                             EOS_Parameter);
        }
    }
    //HLL
    else if (ID_FluxPhaseField == 1){
        FluxPhasePhaseField_HLL(FluxPhase,
                                cmax,
                                ConservativePhaseL,
                                ConservativePhaseR,
                                Grad_Q,
                                EOS_Parameter,
                                ID_BE);
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
#endif