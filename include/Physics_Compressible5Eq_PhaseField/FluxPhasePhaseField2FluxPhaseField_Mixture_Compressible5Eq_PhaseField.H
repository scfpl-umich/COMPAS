// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhasePhaseField2FluxPhaseField_Mixture_Compressible5Eq_PhaseField_H_
#define FluxPhasePhaseField2FluxPhaseField_Mixture_Compressible5Eq_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhasePhaseField2FluxPhaseField_Mixture(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                                                amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*2> const& FluxPhase,
                                                amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeL,
                                                amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeR,
                                                Parm const& EOS_Parameter)
{
    amrex::Real J1_rho1;
    amrex::Real J2_rho2;
    amrex::Real J1_rho1_u1;
    amrex::Real J2_rho2_u2;
#if (AMREX_SPACEDIM > 1)
    amrex::Real J1_rho1_v1;
    amrex::Real J2_rho2_v2;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real J1_rho1_w1;
    amrex::Real J2_rho2_w2;
#endif
    amrex::Real J1_rho1_E1;
    amrex::Real J2_rho2_E2;
    amrex::Real J1;
    amrex::Real J2;
    VectorPhase2Component(J1_rho1, J2_rho2,
                          J1_rho1_u1, J2_rho2_u2,
#if (AMREX_SPACEDIM > 1)
                          J1_rho1_v1, J2_rho2_v2,
#endif
#if (AMREX_SPACEDIM > 2)
                          J1_rho1_w1, J2_rho2_w2,
#endif
                          J1_rho1_E1, J2_rho2_E2,
                          J1, J2,
                          FluxPhase);
    
    FluxPhaseField_Mixture(Flux,
                           ConservativeL,
                           ConservativeR,
                           J1_rho1,
                           J2_rho2,
                           J1,
                           J2,
                           EOS_Parameter);
    
}

#endif