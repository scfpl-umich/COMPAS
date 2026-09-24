// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2Pressure_Compressible6Eq_IE_NPhase_H_
#define Conservative2Pressure_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2Pressure(amrex::Real& Pressure,
                           amrex::Array<amrex::Real,NSTATE> const& Conservative,
                           Parm const& EOS_Parameter)
{
   amrex::Array<amrex::Real, NSTATE> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);

    amrex::Array<amrex::Real,NPHASE> Alpha_rho;
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e;
    AMREX_D_TERM(amrex::Real u;,
                 amrex::Real v;,
                 amrex::Real w);
    amrex::Real p;
    amrex::Array<amrex::Real,NPHASE> Alpha;
    Primitive2State(Alpha,
                    Alpha_rho,
                    Alpha_rho_e,
       AMREX_D_DECL(u, v, w),
                    p,
                    Primitive);
    Pressure = p;
    
}

#endif
#endif