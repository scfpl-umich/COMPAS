// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Viscosity_Face_D_Compressible5Eq_NPhase_H_
#define Conservative2Viscosity_Face_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Viscosity_Face_D(amrex::Real& mu_Face,
                                   amrex::Real& mub_Face,
                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_m,
                                   amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                   Parm const& EOS_Parameter)
{
    amrex::Real mu_m;
    amrex::Real mub_m;
    Conservative2Viscosity(mu_m,
                           mub_m,
                           Conservative_m,
                           EOS_Parameter);

    amrex::Real mu_c;
    amrex::Real mub_c;
    Conservative2Viscosity(mu_c,
                           mub_c,
                           Conservative_c,
                           EOS_Parameter);

    Interpolation2Face_Linear_D(mu_Face,
                                mu_m,
                                mu_c);

    Interpolation2Face_Linear_D(mub_Face,
                                mub_m,
                                mub_c);
    
}

#endif
#endif