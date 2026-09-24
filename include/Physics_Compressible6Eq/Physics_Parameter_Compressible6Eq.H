// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Parameter_Compressible6Eq_H_
#define Physics_Parameter_Compressible6Eq_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>
#if (PHYSICS==SIXEQS)

enum class Riemann : int {
    LLF,
    HLL,
    HLLC
};

struct Physics_Parameter
{
    int relaxation = 0;
    int pressure_relaxation = 0;
    int pressure_temperature_relaxation = 0;
    int pressure_temperature_gibbs_relaxation = 0;
    amrex::Real pressure_relaxation_implicit = 1.0;

    Riemann RiemannSolver = Riemann::HLLC;
    int ID_Riemann = 2;

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("Physics");
        pp.query("pressure_relaxation",pressure_relaxation);
        pp.query("pressure_temperature_relaxation",pressure_temperature_relaxation);
        pp.query("pressure_temperature_gibbs_relaxation",pressure_temperature_gibbs_relaxation);
        pp.query("pressure_relaxation_implicit",pressure_relaxation_implicit);

        if (pressure_temperature_gibbs_relaxation){
            amrex::Abort("Physics.pressure_temperature_gibbs_relaxation is not implemented yet. "
                         "Use Physics.pressure_relaxation or Physics.pressure_temperature_relaxation.");
        }

        relaxation = pressure_relaxation || pressure_temperature_relaxation || pressure_temperature_gibbs_relaxation;
    };
};

#endif
#endif