// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldMechanism_ACAdv_D_H_
#define FDM_PhaseFieldMechanism_ACAdv_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldMechanism_ACAdv_D(amrex::Real& PF_Mechanism,
                                     amrex::Real const& PF_c,
                                     amrex::Real const& PF_xm,
                                     amrex::Real const& PF_xp,
                                     amrex::Real const& Normalx_c,
                                     amrex::Real const& Normalx_xm,
                                     amrex::Real const& Normalx_xp,
#if (AMREX_SPACEDIM > 1)
                                     amrex::Real const& PF_ym,
                                     amrex::Real const& PF_yp,
                                     amrex::Real const& Normaly_c,
                                     amrex::Real const& Normaly_ym,
                                     amrex::Real const& Normaly_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                     amrex::Real const& PF_zm,
                                     amrex::Real const& PF_zp,
                                     amrex::Real const& Normalz_c,
                                     amrex::Real const& Normalz_zm,
                                     amrex::Real const& Normalz_zp,
#endif
                                     amrex::Real const& Mob,
                                     amrex::Real const& Lambda,
                                     amrex::Real const& Eta,
                        AMREX_D_DECL(amrex::Real const& dx,
                                     amrex::Real const& dy,
                                     amrex::Real const& dz),
                                     int ID_Potential = 0)
{
    amrex::Real ChemicalPotential;
    FDM_ChemicalPotential_Adv_D(ChemicalPotential,
                                PF_c,
                                PF_xm, PF_xp,
                                Normalx_c, Normalx_xm, Normalx_xp,
#if (AMREX_SPACEDIM > 1)
                                PF_ym, PF_yp,
                                Normaly_c, Normaly_ym, Normaly_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                PF_zm, PF_zp,
                                Normalz_c, Normalz_zm, Normalz_zp,
#endif
                                Lambda, Eta,
                   AMREX_D_DECL(dx, dy, dz),
                                ID_Potential);
    
    PhaseFieldMechanism_AC(PF_Mechanism,
                           ChemicalPotential, Mob);
    
}

#endif