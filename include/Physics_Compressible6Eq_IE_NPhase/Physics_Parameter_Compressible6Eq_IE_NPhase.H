// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Parameter_Compressible6Eq_IE_NPhase_H_
#define Physics_Parameter_Compressible6Eq_IE_NPhase_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

enum class Riemann : int {
    LLF,
    HLL,
    HLLC
};

struct Physics_Parameter
{
    int pressure_relaxation = 0;
    // Volume fraction at or below which a phase keeps its volume fraction in the pressure
    // relaxation (it reaches the common pressure through its energy); 0 relaxes every present phase
    amrex::Real pressure_relaxation_trace = 0.0;

    Riemann RiemannSolver = Riemann::HLLC;
    int ID_Riemann = 2;

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("Physics");
        pp.query("pressure_relaxation",pressure_relaxation);
        pp.query("pressure_relaxation_trace",pressure_relaxation_trace);
        if (pressure_relaxation_trace < 0.0 || pressure_relaxation_trace >= 1.0){
            amrex::Abort("Physics.pressure_relaxation_trace must be in [0, 1).");
        }
    };
};

#endif
#endif