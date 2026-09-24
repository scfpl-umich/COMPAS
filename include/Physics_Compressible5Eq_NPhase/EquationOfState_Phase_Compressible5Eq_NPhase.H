// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EquationOfState_Phase_Compressible5Eq_NPhase_H_
#define EquationOfState_Phase_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void EquationOfState_Phase(amrex::Real& Output,
                           amrex::Real const& Input,
                           amrex::Real const& Alpha,
                           amrex::Real const& Alpha_rho,
                           Parm const& EOS_Parameter,
                           int ID_Phase,
                           bool PressureBase,
                           bool OutputInternalEnergy)
{
    amrex::Real EOS_AP = EOS_Parameter.A[ID_Phase];
    amrex::Real EOS_BP = EOS_Parameter.B[ID_Phase];
    amrex::Real EOS_D  = EOS_Parameter.D[ID_Phase];

    //Pressure-Based EOS
    if (PressureBase){
        //From Pressure to Alpha_rho_e
        if (OutputInternalEnergy){
            Output = Alpha*EOS_AP*Input + Alpha*EOS_BP + Alpha_rho*EOS_D;
        }
        //From Alpha_rho_e to Alpha_p
        else{
            Output = (Input - Alpha_rho*EOS_D - Alpha*EOS_BP)/EOS_AP;
        }
    }
    //Temperature-Based EOS
    else{
        amrex::Real EOS_CT = EOS_Parameter.cv[ID_Phase];
        amrex::Real EOS_BT = EOS_BP/(1.0 + EOS_AP);
        
        //From Temperature to Alpha_rho_e
        if (OutputInternalEnergy){
            Output = Alpha_rho*EOS_CT*Input + Alpha*EOS_BT + Alpha_rho*EOS_D;
        }
        //From Alpha_rho_e to Alpha_rho_T
        else{
           Output = (Input - Alpha_rho*EOS_D - Alpha*EOS_BT)/EOS_CT;
        }
    }
}

#endif
#endif