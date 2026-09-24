// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.


#ifndef _Adv_K_H_
#define _Adv_K_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <Operators_FV.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_SurfaceIntegral(int i, int j, int k, int v,
                         amrex::Array4<amrex::Real> const& dof_out,
                         amrex::Array4<amrex::Real const> const& dof_in,
            AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                         amrex::Array4<amrex::Real const> const& flxy,
                         amrex::Array4<amrex::Real const> const& flxz),
            AMREX_D_DECL(amrex::Real dx, amrex::Real dy, amrex::Real dz))
{
    amrex::Real Div_f;
    FVM_Divergence_K(i, j, k, v,
                     Div_f,
        AMREX_D_DECL(flxx, flxy, flxz),
        AMREX_D_DECL(dx  , dy  , dz  ));

    dof_out(i,j,k,v) = -Div_f;

}



template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_SurfaceIntegral_Coord(int i, int j, int k,
                            amrex::Array4<amrex::Real> const& dof_out,
                            amrex::Array4<amrex::Real const> const& dof_in,
                            AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                                         amrex::Array4<amrex::Real const> const& flxy,
                                         amrex::Array4<amrex::Real const> const& flxz),
                            int ID_COORDSYS,
                            TDX const & dx,
                            TP const & prob_lo)
{

    amrex::Real tint;
    GetTimeIntegral_K(tint,i,j,k,ID_COORDSYS,dx,prob_lo);
    AMREX_D_TERM(
    amrex::Real fintxm;
    amrex::Real fintxp;
    GetFluxIntegral_K(fintxm,i,j,k,ID_COORDSYS,0,-1.,dx,prob_lo);
    GetFluxIntegral_K(fintxp,i,j,k,ID_COORDSYS,0,1.,dx,prob_lo);,
    amrex::Real fintym;
    amrex::Real fintyp;
    GetFluxIntegral_K(fintym,i,j,k,ID_COORDSYS,1,-1.,dx,prob_lo);
    GetFluxIntegral_K(fintyp,i,j,k,ID_COORDSYS,1,1.,dx,prob_lo);,
    amrex::Real fintzm;
    amrex::Real fintzp;
    GetFluxIntegral_K(fintzm,i,j,k,ID_COORDSYS,2,-1.,dx,prob_lo);
    GetFluxIntegral_K(fintzp,i,j,k,ID_COORDSYS,2,1.,dx,prob_lo););

    for (int v = 0; v < NSTATE; v++){
        dof_out(i,j,k,v) = 
                    ( AMREX_D_TERM( (flxx(i,j,k,v)*fintxm - flxx(i+1,j,k,v)*fintxp),
                                  + (flxy(i,j,k,v)*fintym - flxy(i,j+1,k,v)*fintyp),
                                  + (flxz(i,j,k,v)*fintzm - flxz(i,j,k+1,v)*fintzp)) )/tint;
    }

}




#endif
