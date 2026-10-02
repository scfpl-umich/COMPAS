// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EquationOfState_Mixture_Compressible5Eq_H_
#define EquationOfState_Mixture_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void EquationOfState_Mixture(amrex::Real& Output,
                             amrex::Real const& Input,
                             amrex::Real const& Alpha1,
                             amrex::Real const& Alpha1_rho1,
                             amrex::Real const& Alpha2_rho2,
                             Parm const& EOS_Parameter,
                             bool PressureBase,
                             bool OutputInternalEnergy)
{
    amrex::Real Alpha2 = 1.0 - Alpha1;
    // NASG: phi_k = alpha_k - b_k alpha_k rho_k replaces alpha_k where it multiplies A_k or B_k
    // (phi_k = alpha_k exactly for the other EOS, b_k = 0)
    amrex::Real Phi1 = NASG_Phi(Alpha1, Alpha1_rho1, EOS_Parameter.C1);
    amrex::Real Phi2 = NASG_Phi(Alpha2, Alpha2_rho2, EOS_Parameter.C2);

    amrex::Real EOS_D1 = EOS_Parameter.D1;
    amrex::Real EOS_D2 = EOS_Parameter.D2;
    amrex::Real rho_D = Alpha1_rho1*EOS_D1
                       +Alpha2_rho2*EOS_D2;
        
    //Pressure-Based EOS
    if (PressureBase){
        amrex::Real EOS_AP1 = EOS_Parameter.A1;
        amrex::Real EOS_AP2 = EOS_Parameter.A2;
        amrex::Real AP = Phi1*EOS_AP1
                        +Phi2*EOS_AP2;
        
        amrex::Real EOS_BP1 = EOS_Parameter.B1;
        amrex::Real EOS_BP2 = EOS_Parameter.B2;
        amrex::Real BP = Phi1*EOS_BP1
                        +Phi2*EOS_BP2;
        
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
        amrex::Real EOS_CT1 = EOS_Parameter.cv_1;
        amrex::Real EOS_CT2 = EOS_Parameter.cv_2;
        amrex::Real rho_CT = Alpha1_rho1*EOS_CT1
                            +Alpha2_rho2*EOS_CT2;
        
        amrex::Real EOS_BT1 = EOS_Parameter.B1/(1.0 + EOS_Parameter.A1);
        amrex::Real EOS_BT2 = EOS_Parameter.B2/(1.0 + EOS_Parameter.A2);
        amrex::Real BT = Phi1*EOS_BT1
                        +Phi2*EOS_BT2;
        
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