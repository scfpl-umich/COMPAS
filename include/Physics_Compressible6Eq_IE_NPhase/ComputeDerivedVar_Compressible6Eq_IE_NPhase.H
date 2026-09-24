// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef DERIVEDVAR_Compressible6Eq_IE_NPhase_H_
#define DERIVEDVAR_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)
#include<Conservative2Velocity_Compressible6Eq_IE_NPhase.H>
#include<DerivedVelocityGradient.H>

template<typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_derived_var (int i, 
                     int j, 
                     int k,
                     T dvar,
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
    amrex::Real c_avg = 1.0;
// #if (SPACE == DGALERKIN)
//     c_avg = 1./std::pow(std::sqrt(2.),AMREX_SPACEDIM);
// #endif

    // amrex::Print()<<"c_avg "<<c_avg<<"\n";

    amrex::Array<amrex::Real, NSTATE> Conservative;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    for (int v = 0; v < NSTATE; v++){
        Conservative[v] = statein(i,j,k,v)*c_avg;
    }


    Conservative2Primitive(Primitive,Conservative,parm);

    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;

    amrex::Real rho_u = 0.; 
    amrex::Real rho_v = 0.;
    amrex::Real rho_w = 0.;
    amrex::Real rho_E;

    Conservative2State(Alpha,Alpha_rho,Alpha_rho_e, AMREX_D_DECL(rho_u, rho_v, rho_w),rho_E, Conservative);

   
    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;
    
    amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
    amrex::Real  p;

    Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), p, Primitive);

    amrex::Real rho = 0.0;
    amrex::Real Phase = 0.0;
    amrex::Real PhaseSq = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
        Phase += (1.+i_Phase)*Alpha[i_Phase];
        PhaseSq += Alpha[i_Phase]*Alpha[i_Phase];
    }

    amrex::Real Temp;
    Conservative2Temperature(Temp,Conservative, parm);

    amrex::Real c;
    SoundSpeed(c, Alpha, Alpha_rho, Alpha_Ps,parm);

    amrex::Real Mach = std::sqrt(u*u + v*v + w*w)/c;

    amrex::Real vorticity_mag = 0.0;
    amrex::Real q_criterion = 0.0;
    amrex::Real lambda2 = 0.0;
    bool need_gradient = false;
    for (int var = 0; var < num_output_vars; var++){
        need_gradient = need_gradient || (dvar[var] >= 4*NPHASE+12 && dvar[var] <= 4*NPHASE+14);
    }
    if (need_gradient){
        VelocityGradientInvariants(i,j,k,dx,dy,dz,statein,vorticity_mag,q_criterion,lambda2);
    }

    for(int var = 0; var < num_output_vars; var++){
        amrex::Real value = 0.0;

        if ( dvar[var] <=  NPHASE -1){
            // Volume fractions
            value = Alpha[dvar[var]];
            // amrex::Print()<<"alpha "<<dvar[var]<<" "<<value<<"\n";
        } else if (dvar[var] >  NPHASE -1 && dvar[var] <=  2*NPHASE -1){
            // phasic densities
            value = Alpha_rho[dvar[var] - NPHASE];
        } else if (dvar[var] >  2*NPHASE -1 && dvar[var] <=  3*NPHASE-1){
            // phasic energies
            value = Alpha_rho_e[dvar[var] - 2*NPHASE];
        } else if (dvar[var] >  3*NPHASE -1 && dvar[var] <=  4*NPHASE-1){
            // phasic pressures
            value = Alpha_Ps[dvar[var] - 3*NPHASE];
        }
        switch (dvar[var]){
        // "rho_u","rho_v","rho_w","rho_E","rho","u","v","w","P","T","c","Mach","vorticity_mag","q_criterion","lambda2","Phase","PhaseSq"

            case  4*NPHASE: value = rho_u; break;
            case  4*NPHASE + 1: value = rho_v; break;
            case  4*NPHASE + 2: value = rho_w; break;
            case  4*NPHASE + 3: value = rho_E; break;

            case  4*NPHASE + 4: value = rho; break;
            case  4*NPHASE + 5: value = u; break;
            case  4*NPHASE + 6: value = v; break;
            case  4*NPHASE + 7: value = w; break; 
            case  4*NPHASE + 8: value = p; break;
            case  4*NPHASE + 9: value = Temp; break;
            case  4*NPHASE + 10: value = c; break;
            case  4*NPHASE + 11: value = Mach; break;

            case  4*NPHASE + 12: value = vorticity_mag; break;
            case  4*NPHASE + 13: value = q_criterion; break;
            case  4*NPHASE + 14: value = lambda2; break;
            case  4*NPHASE + 15: value = Phase; break;
            case  4*NPHASE + 16: value = PhaseSq; break;
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
