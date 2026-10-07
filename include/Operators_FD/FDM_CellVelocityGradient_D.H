// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_CellVelocityGradient_D_H_
#define FDM_CellVelocityGradient_D_H_

// =======================================
// Velocity gradient at the center of a cell by central differences, from the velocities of
// its 2*AMREX_SPACEDIM face neighbors: Vel_m[a][b] and Vel_p[a][b] are u_b in the cells at
// -e_a and +e_a. Grad_Vel[a][b] = d u_b / d x_a = (Vel_p[a][b] - Vel_m[a][b])/(2 dx_a).
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void FDM_CellVelocityGradient_D(amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM>& Grad_Vel,
                                amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Vel_m,
                                amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> const& Vel_p,
                   AMREX_D_DECL(amrex::Real const& dx,
                                amrex::Real const& dy,
                                amrex::Real const& dz))
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx, dy, dz)};
    for (int a = 0; a < AMREX_SPACEDIM; ++a){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Grad_Vel[a][b] = (Vel_p[a][b] - Vel_m[a][b])/(2.0*dX[a]);
        }
    }
}

#endif
