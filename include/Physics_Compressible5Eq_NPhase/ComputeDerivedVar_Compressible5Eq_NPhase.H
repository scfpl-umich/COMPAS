// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef DERIVEDVAR_Compressible5Eq_NPhase_H_
#define DERIVEDVAR_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)
#include<Conservative2Velocity_Compressible5Eq_NPhase.H>
#include<DerivedVelocityGradient.H>


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_derived_var (int i, 
                     int j, 
                     int k,
                     // amrex::Array<int, NSTATE> dvar,
                     // amrex::Vector<int> dvar,
                     amrex::Array<int, 2*NSTATE> const dvar,
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

    amrex::Array<amrex::Real, NSTATE> Conservative;
    amrex::Array<amrex::Real, NSTATE> Primitive;
    for (int v = 0; v < NSTATE; v++){
        Conservative[v] = statein(i,j,k,v)*c_avg;
    }

    Conservative2Primitive(Primitive,Conservative,parm);

    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Real rho_u = 0.0;
    amrex::Real rho_v = 0.0;
    amrex::Real rho_w = 0.0;
    amrex::Real rhoE;
    amrex::Array<amrex::Real, NPHASE> Alpha;

    Conservative2State(Alpha_rho,
          AMREX_D_DECL(rho_u,
                       rho_v,
                       rho_w),
                       rhoE,
                       Alpha,
                       Conservative);

    amrex::Real u = 0.0;
    amrex::Real v = 0.0;
    amrex::Real w = 0.0;

    amrex::Real P;
    Primitive2State(Alpha_rho,
        AMREX_D_DECL(u,v,w),
                    P,
                    Alpha,
                    Primitive);

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
    SoundSpeed(c, P, Alpha, Alpha_rho, parm, 0);


    amrex::Real Mach = std::sqrt(u*u + v*v + w*w)/c;

    amrex::Real vorticity_mag = 0.0;
    amrex::Real q_criterion = 0.0;
    amrex::Real lambda2 = 0.0;
    bool need_gradient = false;
    for (int var = 0; var < num_output_vars; var++){
        need_gradient = need_gradient || (dvar[var] >= 2*NPHASE+12 && dvar[var] <= 2*NPHASE+14);
    }
    if (need_gradient){
        VelocityGradientInvariants(i,j,k,dx,dy,dz,statein,vorticity_mag,q_criterion,lambda2);
    }

    for (int var = 0; var < num_output_vars; var++){

        amrex::Real value = 0.0; 
        if (dvar[var] < NPHASE){
            value = Alpha_rho[dvar[var]];
        } else if (dvar[var] > NPHASE + 3 && dvar[var] <= 2*NPHASE + 3){
            value = Alpha[dvar[var] - NPHASE - 4];
        }

        switch (dvar[var]){
        //     // Fill Alpha_rho above
            case NPHASE: value = rho_u; break;
            case NPHASE+1: value = rho_v; break;
            case NPHASE+2: value = rho_w; break;
            case NPHASE+3: value = rhoE; break;
        //     // Fill Alpha above
            case 2*NPHASE+4: value = rho; break;
            case 2*NPHASE+5: value = u; break;
            case 2*NPHASE+6: value = v; break;
            case 2*NPHASE+7: value = w; break;
            case 2*NPHASE+8: value = P; break;
            case 2*NPHASE+9: value = Temp; break;
            case 2*NPHASE+10: value = c; break;
            case 2*NPHASE+11: value = Mach; break;
            case 2*NPHASE+12: value = vorticity_mag; break;
            case 2*NPHASE+13: value = q_criterion; break;  
            case 2*NPHASE+14: value = lambda2; break;
            case 2*NPHASE+15: value = Phase; break;   
            case 2*NPHASE+16: value = PhaseSq; break;
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
