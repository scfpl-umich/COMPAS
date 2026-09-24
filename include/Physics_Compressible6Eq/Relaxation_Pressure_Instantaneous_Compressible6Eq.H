// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Relaxation_Pressure_Instantaneous_Compressible6Eq_H_
#define Relaxation_Pressure_Instantaneous_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Pressure_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
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

    amrex::Real A2 = EOS_Parameter.A2;
    amrex::Real B2 = EOS_Parameter.B2;

    amrex::Real aP = A1*A2 + A1*Alpha1 + A2*Alpha2;
    amrex::Real bP = A1*B2*Alpha1 + A2*B1*Alpha2 - (A2 + 1.0)*A1*Alpha1_p1 - (A1 + 1.0)*A2*Alpha2_p2;
    amrex::Real cP = -B2*A1*Alpha1_p1 - B1*A2*Alpha2_p2;
    amrex::Real pR = (-bP + std::sqrt(bP*bP - 4.0*aP*cP))/(2.0*aP);
    amrex::Real Alpha1R = (A1*Alpha1_p1 + B1*Alpha1 + pR*Alpha1)/((A1 + 1.0)*pR + B1);
    amrex::Real Alpha2R = (A2*Alpha2_p2 + B2*Alpha2 + pR*Alpha2)/((A2 + 1.0)*pR + B2);
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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Pressure_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                       amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                       Parm const& EOS_Parameter,
                                       amrex::Real const& ID_Implicit)
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

    amrex::Real A2 = EOS_Parameter.A2;
    amrex::Real B2 = EOS_Parameter.B2;

    amrex::Real pI = Alpha1_p1 + Alpha2_p2;
    amrex::Real aP = A1*A2 + ID_Implicit*A1*Alpha1 + ID_Implicit*A2*Alpha2;
    amrex::Real bP = A1*(B2 + (1.0 - ID_Implicit)*pI)*Alpha1 + A2*(B1 + (1.0 - ID_Implicit)*pI)*Alpha2 - (A2 + ID_Implicit)*A1*Alpha1_p1 - (A1 + ID_Implicit)*A2*Alpha2_p2;
    amrex::Real cP = -(B2 + (1.0 - ID_Implicit)*pI)*A1*Alpha1_p1 - (B1 + (1.0 - ID_Implicit)*pI)*A2*Alpha2_p2;
    amrex::Real pR = (-bP + std::sqrt(bP*bP - 4.0*aP*cP))/(2.0*aP);
    amrex::Real pIR = (1.0 - ID_Implicit)*pI + ID_Implicit*pR;
    amrex::Real Alpha1R = (A1*Alpha1_p1 + B1*Alpha1 + pIR*Alpha1)/(A1*pR + B1 + pIR);
    amrex::Real Alpha2R = (A2*Alpha2_p2 + B2*Alpha2 + pIR*Alpha2)/(A2*pR + B2 + pIR);
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