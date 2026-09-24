// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef VectorPhase2Vector_Compressible5Eq_NPhase_H_
#define VectorPhase2Vector_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void VectorPhase2Vector(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Vector,
                        amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE> const& VectorPhase)
{
    amrex::Array<amrex::Real,NPHASE> Mass;
    amrex::Array<amrex::Real,NPHASE> Momentum_x;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real,NPHASE> Momentum_y;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real,NPHASE> Momentum_z;
#endif
    amrex::Array<amrex::Real,NPHASE> Energy;
    amrex::Array<amrex::Real,NPHASE> VolumeFraction;
    VectorPhase2Component(Mass,
                          Momentum_x,
#if (AMREX_SPACEDIM > 1)
                          Momentum_y,
#endif
#if (AMREX_SPACEDIM > 2)
                          Momentum_z,
#endif
                          Energy,
                          VolumeFraction,
                          VectorPhase);
    
    amrex::Real TotalMomentum_x = 0.0;
#if (AMREX_SPACEDIM > 1)
    amrex::Real TotalMomentum_y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real TotalMomentum_z = 0.0;
#endif
    amrex::Real TotalEnergy = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        TotalMomentum_x += Momentum_x[i_Phase];
#if (AMREX_SPACEDIM > 1)
        TotalMomentum_y += Momentum_y[i_Phase];
#endif
#if (AMREX_SPACEDIM > 2)
        TotalMomentum_z += Momentum_z[i_Phase];
#endif
        TotalEnergy += Energy[i_Phase];
    }
    Component2Vector(Vector,
                     Mass,
        AMREX_D_DECL(TotalMomentum_x, TotalMomentum_y, TotalMomentum_z),
                     TotalEnergy,
                     VolumeFraction);
    
}

#endif
#endif