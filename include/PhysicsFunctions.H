// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.


#ifndef PHYSICSFUNCTIONS_
#define PHYSICSFUNCTIONS_

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
StiffenedShockConditions(amrex::Real Ma, 
                         amrex::Real u0, 
                         amrex::Real p0, 
                         amrex::Real rho0, 
                         amrex::Real n, 
                         amrex::Real B,
                         T * ret_param)
{
    amrex::Real rho_ps, u_ps, p_ps, rhoE_ps, s_ps;
    amrex::Real c0 = std::sqrt(n*(p0 + B)/rho0);
    p_ps = p0*(1. + 2.*n/(n+1.)*(Ma*Ma - 1.)*(1. + B/p0));
    rho_ps = rho0*((n+1.)*(p_ps+B)/((n-1.)*(p0+B)) + 1.)/((n+1.)/(n-1.) + (p_ps+B)/(p0+B));
    u_ps = u0 + c0/n*(p_ps/p0 - 1.)*p0/(p0 + B)/std::sqrt((n+1.)/(2.*n)*(p_ps/p0 - 1.)*p0/(p0 + B) + 1.);
    s_ps = u0 + c0*std::sqrt((n+1.)/(2.*n)*(p_ps/p0 - 1.)*p0/(p0 + B) + 1.);
    rhoE_ps = rho_ps*u_ps*u_ps/2. + p_ps/(n-1.) + n*B/(n-1.);

    ret_param[0] = rho_ps;
    ret_param[1] = u_ps;
    ret_param[2] = p_ps;
    ret_param[3] = rhoE_ps;
    ret_param[4] = s_ps;

}

#endif /* PHYSICSFUNCTIONS_ */