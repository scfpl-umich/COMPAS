// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef AnalyticNonConservativeFlux_Compressible5Eq_NPhase_H_
#define AnalyticNonConservativeFlux_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)
#include<LoopMacros.H>

template<typename T, typename TS>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AnalyticNonConservativeHyperbolicFlux (int i, 
										int j,
										int k,
										AMREX_D_DECL(
									    T & fxqnc,
									    T & fyqnc,
									    T & fzqnc),
									    TS const & statein,
										T const & uq,
										int source,
										Parm const & EOS_Parameter)
{

	fxqnc = {D_LOOP(NPHASE,0.0),AMREX_D_DECL(0.0,0.0,0.0),0.0,D_LOOP(NPHASE,0.0)};
	fyqnc = {D_LOOP(NPHASE,0.0),AMREX_D_DECL(0.0,0.0,0.0),0.0,D_LOOP(NPHASE,0.0)};
#if (AMREX_SPACEDIM > 2)
	fzqnc = {D_LOOP(NPHASE,0.0),AMREX_D_DECL(0.0,0.0,0.0),0.0,D_LOOP(NPHASE,0.0)};
#endif

    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Conservative;
    amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> Primitive;

    for (int i_State = 0; i_State < 2*NPHASE+AMREX_SPACEDIM+1; i_State++){
        Conservative[i_State] = uq[i_State];
    }

    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);

    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    AMREX_D_TERM( amrex::Real U;,
                  amrex::Real V;,
                  amrex::Real W );
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Primitive2State(Alpha_rho,
       AMREX_D_DECL(U, V, W),
                    p,
                    Alpha,
                    Primitive);

    amrex::Array<amrex::Real,NPHASE> K;                
    NonConservative(K,
                    p,
                    Alpha,
                    EOS_Parameter,
                    source);

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){

        Alpha[i_Phase] = statein(i,j,k,INDEX_VolumeFraction1+i_Phase)/std::pow(std::sqrt(2.),AMREX_SPACEDIM);

    	AMREX_D_TERM(
    	fxqnc[INDEX_VolumeFraction1 + i_Phase] = -Alpha[i_Phase]*(1. + K[i_Phase])*U;,
    	fyqnc[INDEX_VolumeFraction1 + i_Phase] = -Alpha[i_Phase]*(1. + K[i_Phase])*V;,
    	fzqnc[INDEX_VolumeFraction1 + i_Phase] = -Alpha[i_Phase]*(1. + K[i_Phase])*W;);
    }

}

#endif
#endif