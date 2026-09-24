// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef AnalyticHyperbolicFlux_Compressible5Eq_H_
#define AnalyticHyperbolicFlux_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

template<typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
AnalyticHyperbolicFlux (AMREX_D_DECL(
					    T & fxq,
					    T & fyq,
					    T & fzq),
						T const & uq,
						Parm const & parm)
{

	AMREX_D_TERM(
	Conservative2Flux(fxq,uq,parm,0);,
	Conservative2Flux(fyq,uq,parm,1);,
	Conservative2Flux(fzq,uq,parm,2););
	

}

#endif
#endif