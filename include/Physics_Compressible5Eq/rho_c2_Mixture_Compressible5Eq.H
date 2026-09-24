// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef rho_c2_Mixture_Compressible5Eq_H_
#define rho_c2_Mixture_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void rho_c2_Mixture(amrex::Real& rho_c2,
                    amrex::Real const& rho1_c12,
                    amrex::Real const& rho2_c22,
                    amrex::Real const& Alpha1,
                    Parm const& EOS_Parameter,
                    int ID_K)
{
    amrex::Real Alpha2 = 1.0 - Alpha1;
    
    //K_p = 0
    if (ID_K == 0){
        amrex::Real EOS_AP1 = EOS_Parameter.A1;
        amrex::Real EOS_AP2 = EOS_Parameter.A2;
        
        amrex::Real AP_rho_c2 = Alpha1*EOS_AP1*rho1_c12
                               +Alpha2*EOS_AP2*rho2_c22;
        
        amrex::Real AP = Alpha1*EOS_AP1
                        +Alpha2*EOS_AP2;
        
        rho_c2 = AP_rho_c2/AP;
    }
    //K_p = (rho_c2/rhop_cp2-1)
    else if (ID_K == 1){
        rho_c2 = 1.0/( Alpha1/rho1_c12
                      +Alpha2/rho2_c22 );
    }
}

#endif
#endif