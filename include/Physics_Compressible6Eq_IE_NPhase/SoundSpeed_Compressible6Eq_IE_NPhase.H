// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SoundSpeed_Compressible6Eq_IE_NPhase_H_
#define SoundSpeed_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// Frozen sound speed of the mixture, rho c^2 = sum_k alpha_k rho_k c_k^2 (rho_c2_Mixture), with the
// phase relation of rho_c2_Phase for every equation of state (ideal and stiffened gas,
// Mie-Gruneisen, NASG): alpha_k rho_k c_k^2 = ((A_k + 1) alpha_k p_k + alpha_k B_k)/(A_k xi_k),
// with the NASG covolume factor xi_k = 1 - b_k rho_k (exactly 1 when b_k = 0).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
SoundSpeed(amrex::Real& c,
                amrex::Array<amrex::Real,NPHASE> const& Alpha,
                amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
                amrex::Array<amrex::Real,NPHASE> const& Alpha_Ps,
                Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alphap_rhop_cp2;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho_c2_Phase(Alphap_rhop_cp2[i_Phase],
                     Alpha_Ps[i_Phase],
                     Alpha[i_Phase],
                     EOS_Parameter.A[i_Phase],
                     EOS_Parameter.B[i_Phase],
                     NASG_Xi(Alpha[i_Phase], Alpha_rho[i_Phase], EOS_Parameter.C[i_Phase]));
    }
    
    amrex::Real rho_c2;
    rho_c2_Mixture(rho_c2,
                   Alphap_rhop_cp2);
    
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    c = std::sqrt( rho_c2/rho );
}

#endif
#endif
