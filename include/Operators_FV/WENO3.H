// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef WENO3_H_
#define WENO3_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void WENO3(amrex::Array<amrex::Real,2>& f_WENO,
           amrex::Array<amrex::Real,2>& Weight_WENO,
           amrex::Real const& f_m,
           amrex::Real const& f_c,
           amrex::Real const& f_p,
           int ID_Face)
{
    amrex::Real d0;
    amrex::Real d1;
    //f_fp
    if (ID_Face == 1){
        f_WENO[0] = f_c + 0.5*(f_c - f_m);
        f_WENO[1] = f_c + 0.5*(f_p - f_c);
        
        d0 = 1.0/3.0;
    }
    //f_fm
    else if (ID_Face == -1){
        f_WENO[0] = f_c - 0.5*(f_c - f_m);
        f_WENO[1] = f_c - 0.5*(f_p - f_c);
        
        d0 = 2.0/3.0;
    }
    d1 = 1.0 - d0;
    
    amrex::Real IS0 = (f_c - f_m)*(f_c - f_m);
    amrex::Real IS1 = (f_p - f_c)*(f_p - f_c);
    
    amrex::Real IS01;
    amrex::Real IS10;    
    IS01 = IS0/IS1;
    if ( std::isnan(IS01) ){
        IS01 = 1.0;
    }
    IS10 = 1.0/IS01;
    
    Weight_WENO[0] = d0/( d0 + (IS01*IS01)*d1 );
    Weight_WENO[1] = d1/( (IS10*IS10)*d0 + d1 );
}

#endif