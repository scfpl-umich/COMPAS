// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef NonConservative_Compressible5Eq_NPhase_H_
#define NonConservative_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void NonConservative(amrex::Array<amrex::Real,NPHASE>& K,
                     amrex::Real const& p,
                     amrex::Array<amrex::Real,NPHASE> const& Alpha,
                     Parm const& EOS_Parameter,
                     int ID_K)
{
    //K_p = 0
    if (ID_K == 0){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            K[i_Phase] = 0.0;
        }
    }
    //K_p = (rho_c2/rhop_cp2-1)
    else if (ID_K == 1){
        amrex::Array<amrex::Real,NPHASE> rhop_cp2;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            amrex::Real rho_c2_phase;
            rho_c2_Phase(rho_c2_phase,
                         p,
                         EOS_Parameter.A[i_Phase],
                         EOS_Parameter.B[i_Phase]);
            rhop_cp2[i_Phase] = rho_c2_phase;
        }
        
        amrex::Real rho_c2;
        rho_c2_Mixture(rho_c2,
                       rhop_cp2,
                       Alpha,
                       EOS_Parameter,
                       ID_K);
        
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            K[i_Phase] = ( rho_c2/rhop_cp2[i_Phase] - 1.0 );
        }
    }
}

#endif
#endif