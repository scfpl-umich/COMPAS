// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Relaxation_PressureTemperature_Instantaneous_Compressible6Eq_H_
#define Relaxation_PressureTemperature_Instantaneous_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_PressureTemperature_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                                  amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                                  Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real,5+AMREX_SPACEDIM> Primitive;
    Conservative2Primitive(Primitive,
                           Conservative,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1;
    amrex::Real Alpha2_rho2;
    AMREX_D_TERM(amrex::Real u;,
                 amrex::Real v;,
                 amrex::Real w );
    amrex::Real Alpha1_p1;
    amrex::Real Alpha2_p2;
    amrex::Real Alpha1;
    amrex::Real Alpha2;
    Primitive2State(Alpha1_rho1, Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    Alpha1_p1, Alpha2_p2,
                    Alpha1,
                    Primitive);
    Alpha2 = 1.0 - Alpha1;

    amrex::Real A1 = EOS_Parameter.A1;
    amrex::Real B1 = EOS_Parameter.B1;
    amrex::Real C1 = EOS_Parameter.cv_1;
    amrex::Real D1 = EOS_Parameter.D1;

    amrex::Real A2 = EOS_Parameter.A2;
    amrex::Real B2 = EOS_Parameter.B2;
    amrex::Real C2 = EOS_Parameter.cv_2;
    amrex::Real D2 = EOS_Parameter.D2;

    amrex::Real rho_e = A1*Alpha1_p1 + B1*Alpha1 + D1*Alpha1_rho1
                      + A2*Alpha2_p2 + B2*Alpha2 + D2*Alpha2_rho2;
    amrex::Real aPT = C1*Alpha1_rho1 + C2*Alpha2_rho2;
    amrex::Real bPT = C1*D1/A1*Alpha1_rho1*Alpha1_rho1
                    + C2*D2/A2*Alpha2_rho2*Alpha2_rho2
                    + C1*(B1/A1 + B2/(A2 + 1.0))*Alpha1_rho1
                    + C2*(B2/A2 + B1/(A1 + 1.0))*Alpha2_rho2
                    + (C1*D2/A1 + C2*D1/A2)*Alpha1_rho1*Alpha2_rho2
                    - rho_e*(C1/A1*Alpha1_rho1 + C2/A2*Alpha2_rho2);
    amrex::Real cPT = C1*D1/A1*B2/(A2 + 1.0)*Alpha1_rho1*Alpha1_rho1
                    + C2*D2/A2*B1/(A1 + 1.0)*Alpha2_rho2*Alpha2_rho2
                    + (C1*(A1 + 1.0)/A1*Alpha1_rho1 + C2*(A2 + 1.0)/A2*Alpha2_rho2)*B2/(A2 + 1.0)*B1/(A1 + 1.0)
                    + (C1*D2/A1*B2/(A2 + 1.0) + C2*D1/A2*B1/(A1 + 1.0))*Alpha1_rho1*Alpha2_rho2
                    - rho_e*(C1/A1*B2/(A2 + 1.0)*Alpha1_rho1 + C2/A2*B1/(A1 + 1.0)*Alpha2_rho2);
    amrex::Real pR = (-bPT + std::sqrt(bPT*bPT - 4.0*aPT*cPT))/(2.0*aPT);
    amrex::Real Alpha1R = C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0))/(C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0)) + C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0)));
    amrex::Real Alpha2R = C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0))/(C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0)) + C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0)));
    if (std::isnan(Alpha1R) && !std::isnan(Alpha2R)){
        Alpha1R = 1.0 - Alpha2R;
    }
    else if (!std::isnan(Alpha1R) && std::isnan(Alpha2R)){
        Alpha2R = 1.0 - Alpha1R;
    }
    amrex::Real Sum_AlphaR = Alpha1R + Alpha2R;
    Alpha1R /= Sum_AlphaR;
    Alpha2R /= Sum_AlphaR;
    pR = (A1*Alpha1_p1 + A2*Alpha2_p2 + B1*(Alpha1 - Alpha1R) + B2*(Alpha2 - Alpha2R))/(A1*Alpha1R + A2*Alpha2R);
    amrex::Real Alpha1_p1R = Alpha1R*pR;
    amrex::Real Alpha2_p2R = Alpha2R*pR;

    if (std::isnan(Alpha1R) && std::isnan(Alpha2R)){
        Alpha1R    = Alpha1;
        Alpha2R    = Alpha2;
        Alpha1_p1R = Alpha1_p1;
        Alpha2_p2R = Alpha2_p2;
    }

    State2Primitive(Primitive,
                    Alpha1_rho1, Alpha2_rho2,
       AMREX_D_DECL(u, v, w),
                    Alpha1_p1R, Alpha2_p2R,
                    Alpha1R);
    Primitive2Conservative(ConservativeR,
                           Primitive,
                           EOS_Parameter);
    
}

#endif
#endif