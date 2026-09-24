// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef AnalyticNonConservativeFlux_Compressible5Eq_H_
#define AnalyticNonConservativeFlux_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

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
										Parm const & parm)
{
	if (source){

	} else {

		amrex::Real K = -statein(i,j,k,INDEX_VolumeFraction1)/std::pow(std::sqrt(2.),AMREX_SPACEDIM);

		fxqnc = {0.0,0.0,AMREX_D_DECL(0.0,0.0,0.0),0.0,0.0};
		fyqnc = {0.0,0.0,AMREX_D_DECL(0.0,0.0,0.0),0.0,0.0};
#if (AMREX_SPACEDIM > 2)
		fzqnc = {0.0,0.0,AMREX_D_DECL(0.0,0.0,0.0),0.0,0.0};
#endif
		AMREX_D_TERM(
        fxqnc[INDEX_VolumeFraction1] = K*uq[INDEX_MomentumX]/(uq[INDEX_Mass1] + uq[INDEX_Mass2]);,
        fyqnc[INDEX_VolumeFraction1] = K*uq[INDEX_MomentumY]/(uq[INDEX_Mass1] + uq[INDEX_Mass2]);,
        fzqnc[INDEX_VolumeFraction1] = K*uq[INDEX_MomentumZ]/(uq[INDEX_Mass1] + uq[INDEX_Mass2]);
        );
	}
        

}

#endif
#endif