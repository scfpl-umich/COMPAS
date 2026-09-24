// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef DERIVEDVAR_Compressible6Eq_H_
#define DERIVEDVAR_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)
#include<Conservative2Velocity_Compressible6Eq.H>
#include<DerivedVelocityGradient.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_derived_var (int i, 
                     int j, 
                     int k,
                     amrex::Array<int, 2*NSTATE> const dvar,
                     // amrex::Vector<int> dvar,
                     // T & dvar,
                     int num_output_vars,
                     amrex::Real dx,
                     amrex::Real dy,
                     amrex::Real dz,
                     amrex::Array4<amrex::Real> const& varfab,
                     amrex::Array4<amrex::Real const> const& statein,
                     Parm const& parm)
/** Calculates derived variables at a single node given the state, the output
 * array, and the list of variables to be calculated. Need to indicate if 
 * requested derived variable is not available.
 * @param Q Vector to fill with the requested output variables. Should already be sized
 * @param vnames Vector of requested variables
 * @param U Vector of state variables
 */
{
    amrex::Array<amrex::Real, NSTATE> Conservative;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    for (int iState = 0; iState < NSTATE; iState++){
        Conservative[iState] = statein(i,j,k,iState);
    }
    Conservative2Primitive(Primitive,Conservative,parm);

    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;
    amrex::Real Alpha1_p1;
    amrex::Real Alpha2_p2;
    amrex::Real Alpha1;
    amrex::Real T;
    Primitive2State(Alpha1_rho1,Alpha2_rho2,AMREX_D_DECL(u,v,w),Alpha1_p1,Alpha2_p2,Alpha1,Primitive);
    Conservative2Temperature(T,Conservative,parm);
    amrex::Real c; SoundSpeed(c, Alpha1_p1, Alpha2_p2, Alpha1, Alpha1_rho1, Alpha2_rho2, parm);
    amrex::Real Mach = std::sqrt(AMREX_D_TERM(u*u, + v*v, + w*w))/c;

    amrex::Real vorticity_mag = 0.0;
    amrex::Real q_criterion = 0.0;
    amrex::Real lambda2 = 0.0;
    bool need_gradient = false;
    for (int var = 0; var < num_output_vars; var++){
        need_gradient = need_gradient || (dvar[var] >= 18 && dvar[var] <= 20);
    }
    if (need_gradient){
        VelocityGradientInvariants(i,j,k,dx,dy,dz,statein,vorticity_mag,q_criterion,lambda2);
    }

    for (int var = 0; var < num_output_vars; var++){
        amrex::Real value = 0.0;
        switch (dvar[var]){
            //{"a1rho1","a2rho2","rho_u","rho_v","rho_w","a1rho1E1","a2rho2E2","a1","rho","u","v","w","P","a1P1","a2P2","T","c","Mach","vorticity_mag","q_criterion","lambda2"}
            case 0: value = Alpha1_rho1; break;
            case 1: value = Alpha2_rho2; break;
            case 2: value = Conservative[2]; break;
#if (AMREX_SPACEDIM > 1)
            case 3: value = Conservative[3]; break;
#endif
#if (AMREX_SPACEDIM > 2)
            case 4: value = Conservative[4]; break;
#endif
            case 5: value = Conservative[2+AMREX_SPACEDIM]; break;
            case 6: value = Conservative[3+AMREX_SPACEDIM]; break;
            case 7: value = Conservative[4+AMREX_SPACEDIM]; break;

            case 8: value = Alpha1_rho1+Alpha2_rho2; break;
            case 9: value = u; break;
            case 10: value = v; break;
            case 11: value = w; break;
            case 12: value = Alpha1_p1+Alpha2_p2; break;
            case 13: value = Alpha1_p1; break;
            case 14: value = Alpha2_p2; break;
            case 15: value = T; break;
            case 16: value = c; break;
            case 17: value = Mach; break;
            case 18: value = vorticity_mag; break;
            case 19: value = q_criterion; break;
            case 20: value = lambda2; break;
            default:
                value = 0.0;
                break; 
        }
        varfab(i,j,k,var) = value;
    }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_derived_var (int i, 
                     int j, 
                     int k,
                     int dvar,
                     amrex::Real dx,
                     amrex::Real dy,
                     amrex::Real dz,
                     amrex::Array4<amrex::Real> const& varfab,
                     amrex::Array4<amrex::Real const> const& statein,
                     Parm const& parm)
/** Calculates the single derived variable dvar into component 0 of varfab */
{
    amrex::Array<int, 2*NSTATE> dvars;
    dvars[0] = dvar;
    compute_derived_var(i,j,k,dvars,1,dx,dy,dz,varfab,statein,parm);
}

#endif
#endif
