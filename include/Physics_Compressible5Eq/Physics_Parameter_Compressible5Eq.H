// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Parameter_Compressible5Eq_H_
#define Physics_Parameter_Compressible5Eq_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>
#if (SURFACE_TENSION == true)
#include<SurfaceTension_Parameter.H>
#endif
#if (PHYSICS==FIVEEQS)

enum class Riemann : int {
    LLF,
    HLL,
    HLLC
};

struct Physics_Parameter
{
    int         source_term             = 0;
    int         source_term_step_switch = std::numeric_limits<int>::max();
    amrex::Real source_term_time_switch = 0.0;

    Riemann RiemannSolver = Riemann::HLLC;
    int ID_Riemann = 2;

#if (SURFACE_TENSION == true)
    SurfaceTension_Parameter SurfaceTension_Parm;
#endif

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("Physics");
        pp.query("source_term",source_term);
        pp.query("source_term_step_switch",source_term_step_switch);
        pp.query("source_term_time_switch",source_term_time_switch);

#if (SURFACE_TENSION == true)
        SurfaceTension_Parm.Initialize();
#endif
    };
};

#endif
#endif