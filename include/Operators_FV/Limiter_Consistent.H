// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Limiter_Consistent_H_
#define Limiter_Consistent_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Limiter_Consistent(amrex::Array<amrex::Real,NPHASE>& Theta,
                        amrex::Array<amrex::Real,NPHASE> const& Slope)
{
    amrex::Real ep = 1.0e-10;
    
    amrex::Real Sum_Slope = 0.0;
    amrex::Array<amrex::Real,NPHASE> Abs_Slope;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Sum_Slope += Slope[i_Phase];
        Abs_Slope[i_Phase] = std::abs(Slope[i_Phase]);
    }
    int i_Phase_max = std::distance( std::begin(Abs_Slope),std::max_element(std::begin(Abs_Slope),std::end(Abs_Slope)));
    amrex::Real Slope_max = Slope[i_Phase_max];
    
    if ( std::abs(Slope_max)<ep || std::abs(Sum_Slope)<ep ){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Theta[i_Phase] = 1.0;
        }
    }
    else{
        amrex::Real Slope_maxS = Slope_max - Sum_Slope;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Theta[i_Phase] = std::max( std::min(Slope_max/Slope_maxS,1.0),0.0 );
        }
        Theta[i_Phase_max] = std::max( std::min(Slope_maxS/Slope_max,1.0),0.0 );
    }
       
}

#endif