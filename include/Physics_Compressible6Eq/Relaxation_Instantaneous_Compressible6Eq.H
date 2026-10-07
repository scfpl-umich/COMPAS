// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Relaxation_Instantaneous_Compressible6Eq_H_
#define Relaxation_Instantaneous_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

// Iter and Path report the iterations and the path of the relaxation (Relax_Path* in
// RelaxationPressure_Safeguarded.H, plus 10 where a pressure-temperature relaxation fell back to
// the pressure relaxation); the relaxation test case reads them
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                              amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                              Parm const& EOS_Parameter,
                              int& Iter,
                              int& Path)
{
    Iter = 0;
    Path = Relax_PathModel;
    if (EOS_Parameter.Physics_Parm.pressure_temperature_relaxation){
        Relaxation_PressureTemperature_Instantaneous(ConservativeR,
                                                     Conservative,
                                                     EOS_Parameter,
                                                     Iter, Path);
    }
    else if (EOS_Parameter.Physics_Parm.pressure_temperature_gibbs_relaxation){}
    else{
        Relaxation_Pressure_Instantaneous(ConservativeR,
                                          Conservative,
                                          EOS_Parameter,
                                          Iter, Path);
    }
    
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                              amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                              Parm const& EOS_Parameter)
{
    int Iter, Path;
    Relaxation_Instantaneous(ConservativeR, Conservative, EOS_Parameter, Iter, Path);
}

#endif
#endif