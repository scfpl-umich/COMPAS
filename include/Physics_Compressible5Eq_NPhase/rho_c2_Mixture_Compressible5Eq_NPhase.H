// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef rho_c2_Mixture_Compressible5Eq_NPhase_H_
#define rho_c2_Mixture_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void rho_c2_Mixture(amrex::Real& rho_c2,
                    amrex::Array<amrex::Real,NPHASE> const& rhop_cp2,
                    amrex::Array<amrex::Real,NPHASE> const& Alpha,
                    Parm const& EOS_Parameter,
                    int ID_K)
{
    //K_p = 0
    if (ID_K == 0){
        amrex::Real AP_rho_c2 = 0.0;
        amrex::Real AP = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            AP_rho_c2 += Alpha[i_Phase]*EOS_Parameter.A[i_Phase]*rhop_cp2[i_Phase];
            AP        += Alpha[i_Phase]*EOS_Parameter.A[i_Phase];
        }
        
        rho_c2 = AP_rho_c2/AP;
    }
    //K_p = (rho_c2/rhop_cp2-1)
    else if (ID_K == 1){
        amrex::Real rho_c2_Inverse = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            rho_c2_Inverse += Alpha[i_Phase]/rhop_cp2[i_Phase];
        }
        rho_c2 = 1.0/rho_c2_Inverse;
    }
}

#endif
#endif