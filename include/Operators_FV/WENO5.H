// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef WENO5_H_
#define WENO5_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void WENO5(amrex::Array<amrex::Real,3>& f_WENO,
           amrex::Array<amrex::Real,3>& Weight_WENO,
           amrex::Real const& f_mm,
           amrex::Real const& f_m,
           amrex::Real const& f_c,
           amrex::Real const& f_p,
           amrex::Real const& f_pp,
           int ID_Face)
{
    amrex::Real d0;
    amrex::Real d1;
    amrex::Real d2;
    //f_fp
    if (ID_Face == 1){
        f_WENO[0] = f_c + (2.0/3.0)*(f_p - f_c) - (1.0/6.0)*(f_pp - f_p);
        f_WENO[1] = f_c + (1.0/6.0)*(f_c - f_m) + (1.0/3.0)*(f_p - f_c);
        f_WENO[2] = f_c - (1.0/3.0)*(f_m - f_mm) + (5.0/6.0)*(f_c - f_m);
              
        d0 = 3.0/10.0;
        d2 = 1.0/10.0;
    }
    //f_fm
    if (ID_Face == -1){
        f_WENO[0] = f_c - (5.0/6.0)*(f_p - f_c) + (1.0/3.0)*(f_pp - f_p);
        f_WENO[1] = f_c - (1.0/3.0)*(f_c - f_m) - (1.0/6.0)*(f_p - f_c);
        f_WENO[2] = f_c + (1.0/6.0)*(f_m - f_mm) - (2.0/3.0)*(f_c - f_m);
              
        d0 = 1.0/10.0;
        d2 = 3.0/10.0;
    }
    d1 = 1.0 - d0 - d2;
    
    amrex::Real IS0 = (10.0/3.0)*(f_p - f_c)*(f_p - f_c) - (11.0/3.0)*(f_p - f_c)*(f_pp - f_p) + (4.0/3.0)*(f_pp - f_p)*(f_pp - f_p);
    amrex::Real IS1 = -(5.0/3.0)*(f_c - f_m)*(f_p - f_c) + (4.0/3.0)*(f_c - f_m)*(f_c - f_m) + (4.0/3.0)*(f_p - f_c)*(f_p - f_c);
    amrex::Real IS2 = (10.0/3.0)*(f_c - f_m)*(f_c - f_m) - (11.0/3.0)*(f_m - f_mm)*(f_c - f_m) + (4.0/3.0)*(f_m - f_mm)*(f_m - f_mm);
    
    amrex::Real IS01;
    amrex::Real IS10;            
    IS01 = IS0/IS1;
    if ( std::isnan(IS01) ){
        IS01 = 1.0;
    }
    IS10 = 1.0/IS01;
    
    amrex::Real IS02;
    amrex::Real IS20;
    IS02 = IS0/IS2;
    if ( std::isnan(IS02) ){
        IS02 = 1.0;
    }
    IS20 = 1.0/IS02;
    
    amrex::Real IS12;
    amrex::Real IS21;
    IS12 = IS1/IS2;
    if ( std::isnan(IS12) ){
        IS12 = 1.0;
    }
    IS21 = 1.0/IS12;
    
    Weight_WENO[0] = d0/( d0 + (IS01*IS01)*d1 + (IS02*IS02)*d2 );
    Weight_WENO[1] = d1/( (IS10*IS10)*d0 + d1 + (IS12*IS12)*d2 );
    Weight_WENO[2] = d2/( (IS20*IS20)*d0 + (IS21*IS21)*d1 + d2 );
}

#endif