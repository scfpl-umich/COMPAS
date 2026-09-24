// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef ConservativePhase2State_Compressible_TwoPhase_H_
#define ConservativePhase2State_Compressible_TwoPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void ConservativePhase2State(amrex::Real& Alpha1_rho1,
                             amrex::Real& Alpha2_rho2,
                             amrex::Real& Alpha1_rho1_u1,
                             amrex::Real& Alpha2_rho2_u2,
#if (AMREX_SPACEDIM > 1)
                             amrex::Real& Alpha1_rho1_v1,
                             amrex::Real& Alpha2_rho2_v2,
#endif
#if (AMREX_SPACEDIM > 2)
                             amrex::Real& Alpha1_rho1_w1,
                             amrex::Real& Alpha2_rho2_w2,
#endif
                             amrex::Real& Alpha1_rho1_E1,
                             amrex::Real& Alpha2_rho2_E2,
                             amrex::Real& Alpha1,
                             amrex::Real& Alpha2,
                             amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& ConservativePhase)
{
    VectorPhase2Component(Alpha1_rho1,
                          Alpha2_rho2,
                          Alpha1_rho1_u1,
                          Alpha2_rho2_u2,
#if (AMREX_SPACEDIM > 1)
                          Alpha1_rho1_v1,
                          Alpha2_rho2_v2,
#endif
#if (AMREX_SPACEDIM > 2)
                          Alpha1_rho1_w1,
                          Alpha2_rho2_w2,
#endif
                          Alpha1_rho1_E1,
                          Alpha2_rho2_E2,
                          Alpha1,
                          Alpha2,
                          ConservativePhase);
}

#endif
#endif