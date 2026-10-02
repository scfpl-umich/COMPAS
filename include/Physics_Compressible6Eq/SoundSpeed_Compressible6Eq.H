// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SoundSpeed_Compressible6Eq_H_
#define SoundSpeed_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void SoundSpeed(amrex::Real& c,
                amrex::Real const& Alpha1_p1,
                amrex::Real const& Alpha2_p2,
                amrex::Real const& Alpha1,
                amrex::Real const& Alpha1_rho1,
                amrex::Real const& Alpha2_rho2,
                Parm const& EOS_Parameter,
                bool ID_Squared = false)
{
    // NASG covolume factors xi_k = 1 - b_k rho_k (1 when b_k = 0)
    amrex::Real Xi1 = NASG_Xi(Alpha1, Alpha1_rho1, EOS_Parameter.C1);
    amrex::Real Xi2 = NASG_Xi(1.0 - Alpha1, Alpha2_rho2, EOS_Parameter.C2);

    amrex::Real Alpha1_rho1_c12;
    rho_c2_Phase(Alpha1_rho1_c12,
                 Alpha1_p1,
                 Alpha1,
                 EOS_Parameter.A1,
                 EOS_Parameter.B1,
                 Xi1);
    
    amrex::Real Alpha2_rho2_c22;
    rho_c2_Phase(Alpha2_rho2_c22,
                 Alpha2_p2,
                 (1.0 - Alpha1),
                 EOS_Parameter.A2,
                 EOS_Parameter.B2,
                 Xi2);
    
    amrex::Real rho_c2;
    rho_c2_Mixture(rho_c2,
                   Alpha1_rho1_c12,
                   Alpha2_rho2_c22);
    
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