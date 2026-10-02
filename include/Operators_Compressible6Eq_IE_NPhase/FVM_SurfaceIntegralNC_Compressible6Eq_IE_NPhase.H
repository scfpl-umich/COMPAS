// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_SurfaceIntegralNC_Compressible6Eq_IE_NPhase_H_
#define FVM_SurfaceIntegralNC_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS == SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_SurfaceIntegral_NC(int i, int j, int k,
                            amrex::Array4<amrex::Real> const& dof_out,
                            amrex::Array4<amrex::Real const> const& dof_in,
               AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx,
                            amrex::Array4<amrex::Real const> const& flxy,
                            amrex::Array4<amrex::Real const> const& flxz),
               AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxx_NC,
                            amrex::Array4<amrex::Real const> const& flxy_NC,
                            amrex::Array4<amrex::Real const> const& flxz_NC),
               AMREX_D_DECL(amrex::Real dx,
                            amrex::Real dy,
                            amrex::Real dz),
                            Parm const& EOS_Parameter,
                            Physics_Parameter const& Physics_Parm)
{
    #if (AMREX_SPACEDIM > 1)

    amrex::Array<amrex::Real, NSTATE> Conservative;
    amrex::Array<amrex::Real, NSTATE> Primitive;

    for (int vv = 0; vv < NSTATE; ++vv){
        Conservative[vv] = dof_in(i,j,k,vv);
    }

    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;

    amrex::Real rho_u = 0.; 
    amrex::Real rho_v = 0.;
    amrex::Real rho_w = 0.;
    amrex::Real rho_E;

    Conservative2State(Alpha,Alpha_rho,Alpha_rho_e,AMREX_D_DECL(rho_u, rho_v, rho_w),rho_E,Conservative);
    Conservative2Primitive(Primitive, Conservative, EOS_Parameter);

    amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
    amrex::Real u = 0.; 
    amrex::Real v = 0.;
    amrex::Real w = 0.;
    amrex::Real P;

    Primitive2State(Alpha,Alpha_rho,Alpha_Ps,AMREX_D_DECL(u, v, w),P,Primitive);

    // if(i==48 && j==0){
    //     amrex::Print()<<"i "<<i<<" j "<<j<<"\n";
    //     amrex::Print()<<"flxx_NC "<<flxx_NC(i,j,k,0)<<" flxy_NC "<<flxy_NC(i,j,k, 0)<<"\n"; //<<" flxz_NC "<<flxz_NC(i,j,k, 0)<<"\n";
    //     amrex::Print()<<"i +1, j+1, k+1\n";
    //     amrex::Print()<<"flxx_NC "<<flxx_NC(i+1,j,k,0)<<" flxy_NC "<<flxy_NC(i,j+1,k, 0)<<"\n"; //<<" flxz_NC "<<flxz_NC(i,j,k+1, 0)<<"\n";
    // }


    for(int i_phase=0; i_phase<NPHASE; ++i_phase){
        // Alpha
        dof_out(i,j,k, i_phase) -= Alpha[i_phase]*
                ( AMREX_D_TERM( (flxx_NC(i,j,k, 0) - flxx_NC(i+1,j,k, 0)) / dx,
                              + (flxy_NC(i,j,k, 0) - flxy_NC(i,j+1,k, 0)) / dy,
                              + (flxz_NC(i,j,k, 0) - flxz_NC(i,j,k+1, 0)) / dz ) );

        //Alpha_rho_e
        dof_out(i,j,k, 2*NPHASE + i_phase) += Alpha_Ps[i_phase]*
                ( AMREX_D_TERM( (flxx_NC(i,j,k, 0) - flxx_NC(i+1,j,k, 0)) / dx,
                              + (flxy_NC(i,j,k, 0) - flxy_NC(i,j+1,k, 0)) / dy,
                              + (flxz_NC(i,j,k, 0) - flxz_NC(i,j,k+1, 0)) / dz ) );
    
        
    }

    #endif

#if (DIFFUSION == true)
    // Viscous dissipation of each phase internal energy, alpha_k Phi_k with Phi_k = tau_k : grad u
    // (Huang, arXiv:2609.18085, Eq. 20). The conservative viscous flux u.tau goes to rho E
    // (FluxDiffusion_Compressible6Eq_IE_NPhase.H).
    amrex::Array<amrex::Real,NPHASE> Phi_Phase;
    FDM_Conservative2ViscousDissipation_K(i, j, k,
                                          Phi_Phase,
                                          dof_in,
                             AMREX_D_DECL(dx, dy, dz),
                                          EOS_Parameter);
    for (int i_phase = 0; i_phase < NPHASE; ++i_phase){
        dof_out(i,j,k, INDEX_Energy1 + i_phase) += dof_in(i,j,k, INDEX_VolumeFraction1 + i_phase)*Phi_Phase[i_phase];
    }
#endif

}


#endif
#endif