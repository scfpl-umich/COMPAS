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
                    amrex::Real const& Alpha1_rho1,
                    amrex::Real const& Alpha2_rho2,
                    Parm const& EOS_Parameter,
                    int ID_K)
{
    amrex::Real Alpha2 = 1.0 - Alpha1;
    // rho1_c12 and rho2_c22 are gamma_k (p + pinf_k) = (1 - b_k rho_k) rho_k c_k^2 (rho_c2_Phase).
    // NASG: phi_k = alpha_k - b_k alpha_k rho_k (alpha_k when b_k = 0) in the A sum (K = 0) and
    // in Wood's sum, sum alpha_k/(rho_k c_k^2) = sum phi_k/(gamma_k (p + pinf_k))
    amrex::Real Phi1 = NASG_Phi(Alpha1, Alpha1_rho1, EOS_Parameter.C1);
    amrex::Real Phi2 = NASG_Phi(Alpha2, Alpha2_rho2, EOS_Parameter.C2);
    
    //K_p = 0
    if (ID_K == 0){
        amrex::Real EOS_AP1 = EOS_Parameter.A1;
        amrex::Real EOS_AP2 = EOS_Parameter.A2;
        
        amrex::Real AP_rho_c2 = Alpha1*EOS_AP1*rho1_c12
                               +Alpha2*EOS_AP2*rho2_c22;
        
        amrex::Real AP = Phi1*EOS_AP1
                        +Phi2*EOS_AP2;
        
        rho_c2 = AP_rho_c2/AP;
    }
    //K_p = (rho_c2/rhop_cp2-1)
    else if (ID_K == 1){
        rho_c2 = 1.0/( Phi1/rho1_c12
                      +Phi2/rho2_c22 );
    }
}

#endif
#endif