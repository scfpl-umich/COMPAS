// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EquationOfState_Phase_Compressible6Eq_IE_NPhase_H_
#define EquationOfState_Phase_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

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
    amrex::Real EOS_AP;
    amrex::Real EOS_D;
    amrex::Real EOS_CT;
    amrex::Real EOS_BP;
    amrex::Real EOS_BT;

    // NASG: alpha rho e = A xi (alpha p) + phi B + alpha rho q, with xi = 1 - b rho and phi = xi alpha
    // (xi = 1 and phi = alpha exactly when b = 0)
    amrex::Real Xi  = NASG_Xi(Alpha, Alpha_rho, EOS_Parameter.C[ID_Phase]);
    amrex::Real Phi = Xi*Alpha;

    EOS_AP = EOS_Parameter.A[ID_Phase];
    EOS_D = EOS_Parameter.D[ID_Phase];
    EOS_CT = EOS_Parameter.cv[ID_Phase];
    EOS_BP = Phi*EOS_Parameter.B[ID_Phase];
    EOS_BT = Phi*EOS_Parameter.B[ID_Phase]/(1.0 + EOS_Parameter.A[ID_Phase]);

    amrex::Real Alpha_rho_D = Alpha_rho*EOS_D;

    // amrex::Print()<<"A1 "<<EOS_Parameter.A[0]<<" A2 "<<EOS_Parameter.A[1]<<"\n";
    // amrex::Print()<<"B1 "<<EOS_Parameter.B[0]<<" B2 "<<EOS_Parameter.B[1]<<"\n";
    // amrex::Print()<<"D1 "<<EOS_Parameter.D[0]<<" D2 "<<EOS_Parameter.D[1]<<"\n";
    // amrex::Abort();
        
    //Pressure-Based EOS
    if (PressureBase){
        amrex::Real BP = EOS_BP;
        
        //From Pressure to Internal Energy
        if (OutputInternalEnergy){
            Output = EOS_AP*Xi*Input + BP + Alpha_rho_D; // Input is alpha_P

        }
        //From Internal Energy to Pressure
        else {
            
            Output = (Input - Alpha_rho_D - BP)/(EOS_AP*Xi); // Input is alpha_rho_e_h

        }
    }
    //Temperature-Based EOS
    else {
        amrex::Real BT = EOS_BT;
        
        //From Temperature to Internal Energy
        if (OutputInternalEnergy){
            Output = EOS_CT*Input + BT + Alpha_rho_D; // Input is alpha_rho_T
        }
        //From Internal Energy to Temperature
        else {
           Output = (Input - Alpha_rho_D - BT)/EOS_CT; // Input is alpha_rho_e_h
        }
    }
    
}

#endif
#endif