// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EquationOfState_Mixture_Compressible5Eq_NPhase_H_
#define EquationOfState_Mixture_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void EquationOfState_Mixture(amrex::Real& Output,
                             amrex::Real const& Input,
                             amrex::Array<amrex::Real,NPHASE> const& Alpha,
                             amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
                             Parm const& EOS_Parameter,
                             bool PressureBase,
                             bool OutputInternalEnergy)
{
    amrex::Real rho_D = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho_D += Alpha_rho[i_Phase]*EOS_Parameter.D[i_Phase];
    }
        
    //Pressure-Based EOS
    if (PressureBase){
        amrex::Real AP = 0.0;
        amrex::Real BP = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            // NASG: phi_k = alpha_k - b_k alpha_k rho_k (alpha_k when b_k = 0)
            amrex::Real Phi = NASG_Phi(Alpha[i_Phase], Alpha_rho[i_Phase], EOS_Parameter.C[i_Phase]);
            AP += Phi*EOS_Parameter.A[i_Phase];
            BP += Phi*EOS_Parameter.B[i_Phase];
        }
        
        //From Pressure to Internal Energy
        if (OutputInternalEnergy){
            Output = AP*Input + BP + rho_D;
        }
        //From Internal Energy to Pressure
        else{
            Output = (Input - rho_D - BP)/AP;
        }
    }
    //Temperature-Based EOS
    else{
        amrex::Real rho_CT = 0.0;
        amrex::Real BT = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            rho_CT += Alpha_rho[i_Phase]*EOS_Parameter.cv[i_Phase];
            amrex::Real Phi = NASG_Phi(Alpha[i_Phase], Alpha_rho[i_Phase], EOS_Parameter.C[i_Phase]);
            BT += Phi*EOS_Parameter.B[i_Phase]/(1.0 + EOS_Parameter.A[i_Phase]);
        }
        
        //From Temperature to Internal Energy
        if (OutputInternalEnergy){
            Output = rho_CT*Input + BT + rho_D;
        }
        //From Internal Energy to Temperature
        else{
           Output = (Input - rho_D - BT)/rho_CT;
        }
    }
}

#endif
#endif