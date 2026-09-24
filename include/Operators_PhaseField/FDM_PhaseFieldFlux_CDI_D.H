// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_CDI_D_H_
#define FDM_PhaseFieldFlux_CDI_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_CDI_D(amrex::Real& PF_Flux,
                              amrex::Real const& PF_c,
                              amrex::Real const& PF_m,
                              amrex::Real const& Normal_c,
                              amrex::Real const& Normal_m,
                              amrex::Real const& Mob,
                              amrex::Real const& Eta,
                              amrex::Real const& dx)
{
    amrex::Real Grad_PF;
    FDM_Gradient_D(Grad_PF,
                   PF_m, PF_c,
                   dx);

    amrex::Real PF;
    Interpolation2Face_Linear_D(PF,
                                PF_m, PF_c);
    
    amrex::Real Normal;
    Interpolation2Face_Linear_D(Normal,
                                Normal_m, Normal_c);
    
    PhaseFieldFlux_CDI(PF_Flux,
                       Grad_PF, PF, Normal,
                       Mob, Eta);
}

#endif