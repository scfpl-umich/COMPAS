// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef EnergyCorrection_Compressible6Eq_IE_NPhase_H_
#define EnergyCorrection_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void EnergyCorrection(amrex::Array<amrex::Real,NSTATE>& ConservativeC,
                      amrex::Array<amrex::Real,NSTATE> const& Conservative,
                      Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,NPHASE> Alpha;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real u;,
                  amrex::Real v;,
                  amrex::Real w );

    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e;

    AMREX_D_TERM( amrex::Real rho_u;, 
                  amrex::Real rho_v;, 
                  amrex::Real rho_w);        
    amrex::Real rho_E;
    Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u, rho_v, rho_w),rho_E, Conservative);

    // amrex::Print()<<"Before correction: ";
    // for(int ip=0; ip<NSTATE; ++ip){
    //     amrex::Print()<<Conservative[ip]<<" ";
    // }
    // amrex::Print()<<"\n";


    amrex::Real sumAlpha = 0.0;

    // Step 1: Clip negative volume fractions
    for (int ip = 0; ip < NPHASE; ++ip) {

        Alpha[ip] = amrex::max(Alpha[ip], 0.0);
        sumAlpha += Alpha[ip];
    }

    // Step 2: Renormalize so that sum(Alpha) = 1

    if (sumAlpha > 0.0) {

        const amrex::Real invSumAlpha = 1.0 / sumAlpha;
        for (int ip = 0; ip < NPHASE; ++ip) {
            Alpha[ip] *= invSumAlpha;
        }
    }

    // amrex::Print()<<"Corrected vofs : "<<Alpha[0]<<" "<<Alpha[1]<<" "<<Alpha[2]<<"\n";



    bool PressureBase = true;//Pressure
    bool OutputInternalEnergy = true;//From Alpha_p to Alpha_rho_e

    amrex::Real rho_e_min = 0.0;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e_min = {0.};
    amrex::Real rho_e = 0.0;

    for (int ip = 0; ip < NPHASE; ++ip){

        EquationOfState_Phase(Alpha_rho_e_min[ip],
                    Alpha[ip]*EOS_Parameter.p_min[ip],
                    Alpha[ip],
                    Alpha_rho[ip],
                    EOS_Parameter,
                    ip,
                    PressureBase,
                    OutputInternalEnergy);

        rho_e_min += Alpha_rho_e_min[ip];
        rho_e += Alpha_rho_e[ip];

        // amrex::Print()<<" Alpha_rho_e_min "<<Alpha_rho_e_min[ip]<<"\n";

    }

    if (rho_e < rho_e_min){
        amrex::Abort("rho_e < rho_e_min");
    }
    else{

        amrex::Real sum_excess = 0.0;

        // ---------------------------------------------------------

        // Compute total minimum energy and available original excess

        // ---------------------------------------------------------

        for (int ip = 0; ip < NPHASE; ++ip) {

            sum_excess += amrex::max( Alpha_rho_e[ip] - Alpha_rho_e_min[ip], 0.0);

        }

        const amrex::Real available = rho_e - rho_e_min;

        // ---------------------------------------------------------

        // Redistribute energy while respecting every phase minimum

        // ---------------------------------------------------------

        if (sum_excess > 0.0) {

            // amrex::Print() <<"rho_e: "<<rho_e<<" rho_e_min: "<<rho_e_min<<"\n";

            // amrex::Print() <<"available: "<<available<<" sum_excess: "<<sum_excess<<"\n";

            const amrex::Real scale = available / sum_excess;

            for (int ip = 0; ip < NPHASE; ++ip) {

                const amrex::Real excess = amrex::max(Alpha_rho_e[ip] - Alpha_rho_e_min[ip],0.0);

                // amrex::Print()<<"scale "<<scale<<" excess "<<excess<<"\n";

                Alpha_rho_e[ip] = Alpha_rho_e_min[ip]+ scale * excess;

            }
        }

        else {

            // Pathological case:

            // all original phasic energies are at or below their minima.

            //

            // rho_e >= rho_e_min, so distribute the remaining energy

            // equally above the minima.

            const amrex::Real extra = available / static_cast<amrex::Real>(NPHASE);

            for (int ip = 0; ip < NPHASE; ++ip) {
                Alpha_rho_e[ip] = Alpha_rho_e_min[ip] + extra;
            }
        }
    }

    // amrex::Print()<<"Corrected energies: "<<Alpha_rho_e[0]<<" "<<Alpha_rho_e[1]<<" "<<Alpha_rho_e[2]<<"\n";

    State2Conservative(ConservativeC, Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u,rho_v,rho_w), rho_E);

    // amrex::Print()<<"After correction: ";
    // for(int ip=0; ip<NSTATE; ++ip){
    //     amrex::Print()<<ConservativeC[ip]<<" ";
    // }
    // amrex::Print()<<"\n";
    
}

#endif
#endif