// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_PhaseFieldFlux_CDI_CompressionInterpolation_D_H_
#define FDM_PhaseFieldFlux_CDI_CompressionInterpolation_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_PhaseFieldFlux_CDI_CompressionInterpolation_D(amrex::Real& PF_Flux,
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

    amrex::Real Compression_c = PF_c*(1.0 - PF_c)*Normal_c;
    amrex::Real Compression_m = PF_m*(1.0 - PF_m)*Normal_m;
    amrex::Real Compression;
    Interpolation2Face_Linear_D(Compression,
                                Compression_m, Compression_c);
    
    PF_Flux = Mob*(Grad_PF - (std::sqrt(2.0)/Eta)*Compression);
}

#endif