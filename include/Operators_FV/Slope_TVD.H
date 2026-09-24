// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Slope_TVD_H_
#define Slope_TVD_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Slope_TVD(amrex::Real& Slope,
               amrex::Real const& Slope_m,
               amrex::Real const& Slope_c,
               amrex::Real const& Slope_p,
               int ID_Limiter)
{
    //Minmod
    if (ID_Limiter == 1){
        Limiter_Minmod(Slope,
                       Slope_m,
                       Slope_p);
    }
    //vanLeer
    else if (ID_Limiter == 2){
        Limiter_vanLeer(Slope,
                        Slope_m,
                        Slope_p);
    }
    //Superbee
    else if (ID_Limiter == 3){
        Limiter_Superbee(Slope,
                         Slope_m,
                         Slope_p);
    }
    //MC
    else if (ID_Limiter == 4){
        Limiter_MC(Slope,
                   Slope_m,
                   Slope_p,
                   Slope_c);
    }
}

#endif