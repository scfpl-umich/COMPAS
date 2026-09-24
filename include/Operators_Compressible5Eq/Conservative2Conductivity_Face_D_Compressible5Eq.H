// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Conductivity_Face_D_Compressible5Eq_H_
#define Conservative2Conductivity_Face_D_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Conductivity_Face_D(amrex::Real& Kappa_Face,
                                      amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_m,
                                      amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_c,
                                      Parm const& EOS_Parameter)
{
    amrex::Real Kappa_m;
    Conservative2Conductivity(Kappa_m,
                              Conservative_m,
                              EOS_Parameter);

    amrex::Real Kappa_c;
    Conservative2Conductivity(Kappa_c,
                              Conservative_c,
                              EOS_Parameter);

    Interpolation2Face_Linear_D(Kappa_Face,
                                Kappa_m,
                                Kappa_c);
    
}

#endif
#endif