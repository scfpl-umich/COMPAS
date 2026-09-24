// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseFieldFlux2GradQPF_H_
#define PhaseFieldFlux2GradQPF_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseFieldFlux2GradQPF(amrex::Real& Grad_QPF,
                            amrex::Real const& PF_Flux,
                            amrex::Real const& PF,
                            amrex::Real const& epsilon)
{
    amrex::Real Flux_PF;
    if (PF < epsilon || PF > 1.0 - epsilon){
        Flux_PF = 0.0;
    }
    else{
        Flux_PF = PF_Flux;
    }
    
    amrex::Real Weight_PF = PF*(1.0 - PF);
    
    Grad_QPF = Flux_PF/(Weight_PF + epsilon);
}

#endif