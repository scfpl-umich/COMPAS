// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SoundSpeed_Compressible5Eq_H_
#define SoundSpeed_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SoundSpeed(amrex::Real& c,
                amrex::Real const& p,
                amrex::Real const& Alpha1,
                amrex::Real const& Alpha1_rho1,
                amrex::Real const& Alpha2_rho2,
                Parm const& EOS_Parameter,
                int ID_K,
                bool ID_Squared = false)
{
    amrex::Real rho1_c12;
    rho_c2_Phase(rho1_c12,
                 p,
                 EOS_Parameter.A1,
                 EOS_Parameter.B1);
    
    amrex::Real rho2_c22;
    rho_c2_Phase(rho2_c22,
                 p,
                 EOS_Parameter.A2,
                 EOS_Parameter.B2);
    
    amrex::Real rho_c2;
    rho_c2_Mixture(rho_c2,
                   rho1_c12,
                   rho2_c22,
                   Alpha1,
                   EOS_Parameter,
                   ID_K);
    
    amrex::Real rho = Alpha1_rho1
                     +Alpha2_rho2;

    if (ID_Squared){
        c = rho_c2/rho;
    } else {
        c = std::sqrt( rho_c2/rho );    
    }
}

#endif
#endif