// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef WeightFunction_Quadratic_H_
#define WeightFunction_Quadratic_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void WeightFunction_Quadratic(amrex::Real& f,
                              amrex::Real const& x,
                              int ID_Derivative)
{
    //f
    if (ID_Derivative == 0){
        f = x*(1.0 - x);
        // f = 1.0;
    }
    //df_dx
    else if (ID_Derivative == 1){
        f = 1.0 - 2.0*x;
    }
    //d2f_dx2
    else if (ID_Derivative == 2){
        f = -2.0;
    }
    else{
        f = 0.0;
    }
    
}

#endif