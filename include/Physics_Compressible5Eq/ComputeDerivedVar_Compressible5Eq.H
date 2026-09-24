// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef DERIVEDVAR_Compressible5Eq_H_
#define DERIVEDVAR_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)
#include<Conservative2Temperature_Compressible5Eq.H>
#include<Conservative2Velocity_Compressible5Eq.H>
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
    amrex::Real p;
    amrex::Real Alpha1;
    amrex::Real T;
    Primitive2State(Alpha1_rho1,Alpha2_rho2,AMREX_D_DECL(u,v,w),p,Alpha1,Primitive);
    Conservative2Temperature(T,Conservative,parm);
    amrex::Real c; SoundSpeed(c,p,Alpha1,Alpha1_rho1,Alpha2_rho2,parm,0);
    amrex::Real Mach = std::sqrt(AMREX_D_TERM(u*u, + v*v, + w*w))/c;

    amrex::Real vorticity_mag = 0.0;
    amrex::Real q_criterion = 0.0;
    amrex::Real lambda2 = 0.0;
    bool need_gradient = false;
    for (int var = 0; var < num_output_vars; var++){
        need_gradient = need_gradient || (dvar[var] >= 15 && dvar[var] <= 17);
    }
    if (need_gradient){
        VelocityGradientInvariants(i,j,k,dx,dy,dz,statein,vorticity_mag,q_criterion,lambda2);
    }

    for (int var = 0; var < num_output_vars; var++){
        amrex::Real value = 0.0;
        switch (dvar[var]){
            //{"a1rho1","a2rho2","rho_u","rho_v","rho_w","rho_E","a1","rho","u","v","w","P","T","c","Mach","vorticity_mag","q_criterion","lambda2"}
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

            case 7: value = Alpha1_rho1+Alpha2_rho2; break;
            case 8: value = u; break;
            case 9: value = v; break;
            case 10: value = w; break;
            case 11: value = p; break;
            case 12: value = T; break;
            case 13: value = c; break;
            case 14: value = Mach; break;
            case 15: value = vorticity_mag; break;
            case 16: value = q_criterion; break;
            case 17: value = lambda2; break;
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
