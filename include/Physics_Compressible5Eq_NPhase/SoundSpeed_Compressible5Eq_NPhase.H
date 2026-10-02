// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SoundSpeed_Compressible5Eq_NPhase_H_
#define SoundSpeed_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SoundSpeed(amrex::Real& c,
                amrex::Real const& p,
                amrex::Array<amrex::Real,NPHASE> const& Alpha,
                amrex::Array<amrex::Real,NPHASE> const& Alpha_rho,
                Parm const& EOS_Parameter,
                int ID_K,
                int ID_Squared = false)
{
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
                   Alpha_rho,
                   EOS_Parameter,
                   ID_K);
    
    amrex::Real rho = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }
    if (ID_Squared){
        c = rho_c2/rho;
    } else {
        c = std::sqrt( rho_c2/rho );    
    }
    
}

#endif
#endif