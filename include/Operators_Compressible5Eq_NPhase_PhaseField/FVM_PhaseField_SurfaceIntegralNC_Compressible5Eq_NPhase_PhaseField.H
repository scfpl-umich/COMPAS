// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_PhaseField_SurfaceIntegralNC_Compressible5Eq_NPhase_PhaseField_H_
#define FVM_PhaseField_SurfaceIntegralNC_Compressible5Eq_NPhase_PhaseField_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PhaseField_FVM_SurfaceIntegral_NC(int i, int j, int k,
                                    amrex::Array4<amrex::Real> const& dof_out,
                                    amrex::Array4<amrex::Real const> const& dof_in,
                       AMREX_D_DECL(amrex::Array4<amrex::Real const> const& flxxNC,
                                    amrex::Array4<amrex::Real const> const& flxyNC,
                                    amrex::Array4<amrex::Real const> const& flxzNC),
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                                    Parm const & EOS_Parameter)
{
}

#endif
#endif