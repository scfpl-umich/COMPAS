// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PRECONSTRUCTORS_H_
#define PRECONSTRUCTORS_H_

#include <Macros.H>

#if (PHYSICS == FIVEEQS_NPHASE)
#define PRECONSTRUCTOR
AMREX_FORCE_INLINE
void 
PreConstructor(){
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesC.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1));
    }
    AMREX_D_TERM(
    vnamesC.push_back("rho_u");,
    vnamesC.push_back("rho_v");,
    vnamesC.push_back("rho_w"));

    vnamesC.push_back("rho_E");
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesC.push_back("a"+std::to_string(i_Phase+1));
    }



    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesD.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1));
    }
    vnamesD.push_back("rho_u");
    vnamesD.push_back("rho_v");
    vnamesD.push_back("rho_w");
    vnamesD.push_back("rho_E");
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesD.push_back("a"+std::to_string(i_Phase+1));
    }

    vnamesD.push_back("rho");
    vnamesD.push_back("u");
    vnamesD.push_back("v");
    vnamesD.push_back("w");
    vnamesD.push_back("P");
    vnamesD.push_back("T");
    vnamesD.push_back("c");
    vnamesD.push_back("Mach");
    vnamesD.push_back("vorticity_mag");
    vnamesD.push_back("q_criterion");
    vnamesD.push_back("lambda2");
    vnamesD.push_back("Phase");
    vnamesD.push_back("PhaseSq");
}

#endif


#if (PHYSICS == SIXEQS_IE_NPHASE)
#define PRECONSTRUCTOR
AMREX_FORCE_INLINE
void 
PreConstructor(){

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesC.push_back("a"+std::to_string(i_Phase+1));
        vnamesD.push_back("a"+std::to_string(i_Phase+1));
    }
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesC.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1));
        vnamesD.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1));
    }

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesC.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1)+"e"+std::to_string(i_Phase+1));
        vnamesD.push_back("a"+std::to_string(i_Phase+1)+"rho"+std::to_string(i_Phase+1)+"e"+std::to_string(i_Phase+1));
    }

    
    // The phasic pressures are derived, not state components
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        vnamesD.push_back("a"+std::to_string(i_Phase+1)+"p"+std::to_string(i_Phase+1));
    }

    AMREX_D_TERM(
    vnamesC.push_back("rho_u");,
    vnamesC.push_back("rho_v");,
    vnamesC.push_back("rho_w"));

    vnamesD.push_back("rho_u");
    vnamesD.push_back("rho_v");
    vnamesD.push_back("rho_w");

    vnamesC.push_back("rho_E");
    vnamesD.push_back("rho_E");

    vnamesD.push_back("rho");
    vnamesD.push_back("u");
    vnamesD.push_back("v");
    vnamesD.push_back("w");
    vnamesD.push_back("P");
    vnamesD.push_back("T");
    vnamesD.push_back("c");
    vnamesD.push_back("Mach");
    vnamesD.push_back("vorticity_mag");
    vnamesD.push_back("q_criterion");
    vnamesD.push_back("lambda2");
    vnamesD.push_back("Phase");
    vnamesD.push_back("PhaseSq");
}

#endif


#endif