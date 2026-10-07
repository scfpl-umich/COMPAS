// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Relaxation_PressureTemperature_Instantaneous_Compressible6Eq_H_
#define Relaxation_PressureTemperature_Instantaneous_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

// Instantaneous pressure-temperature relaxation, in closed form (Huang 2026, Remark 3.7.1): a
// quadratic for the relaxed pressure, the volume fractions at the common temperature, and the
// pressure again from the conservation of E_1 + E_2. The root of the quadratic a p^2 + b p + c = 0
// is the larger one, computed as -2c/(b + sqrt(b^2 - 4ac)) when b > 0: (-b + sqrt(b^2 - 4ac))/(2a)
// loses the digits of a root that is small against the pinf_k scale (4ac/b^2 is 1e-3 to 1e-5 for a
// gas and a liquid near 1 bar, and 1e-9 for a gas near vacuum, where it gave temperatures off by
// 100 % and, with NASG, a liquid at its covolume limit). The closed form is kept wherever it has a
// pressure-temperature equilibrium: alpha_k rho_k >= 0, a finite result with the pressure above
// -pinf_k and a positive temperature for every phase with mass, and (NASG) every phase with mass
// below its covolume limit, alpha_kR > b_k alpha_k rho_k. There is no check of the volume fractions
// against [0, 1]: where these conditions hold, they lie in [0, 1] (alpha_kR = W_k/(W_1 + W_2) with
// W_k = w_k (P + pinf_l) >= 0, l the other phase, for the stiffened gas). Where the equilibrium does
// not exist (a mixture internal energy too low for a positive common temperature, a negative mass,
// or no free volume left by the covolumes), the cell is relaxed in pressure only
// (Relaxation_Pressure_Instantaneous, Path + 10).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_PressureTemperature_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                                  amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                                  Parm const& EOS_Parameter,
                                                  int& Iter,
                                                  int& Path)
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

    amrex::Real Alpha1R;
    amrex::Real Alpha2R;
    amrex::Real Alpha1_p1R;
    amrex::Real Alpha2_p2R;
    // The relaxed pressure, from the conservation of E_1 + E_2
    amrex::Real P_Relaxed;
    if (EOS_Parameter.C1 > 0.0 || EOS_Parameter.C2 > 0.0){
        // NASG (b_1 or b_2 > 0). At common p and T, v_k = cv_k T/(A_k (p + pinf_k)) + b_k and
        // e_k = cv_k T (A_k p + B_k)/(A_k (p + pinf_k)) + q_k. With w_k = alpha_k rho_k cv_k/A_k,
        // Phi = 1 - b_1 alpha_1 rho_1 - b_2 alpha_2 rho_2 and E = rho e - sum alpha_k rho_k q_k, the volume
        // constraint T sum w_k/(p + pinf_k) = Phi and the energy E = T sum w_k (A_k p + B_k)/(p + pinf_k)
        // give a p^2 + b p + c = 0 with the coefficients below (the stiffened-gas ones at Phi = 1), and
        // alpha_1R = b_1 alpha_1 rho_1 + Phi w_1 (p + pinf_2)/(w_1 (p + pinf_2) + w_2 (p + pinf_1)).
        // (the locals C1 and C2 are cv_1 and cv_2; the covolumes are EOS_Parameter.C1 and C2).
        // E comes from the conservative state, without the covolume floor (see the pressure relaxation).
        amrex::Real b1 = EOS_Parameter.C1;
        amrex::Real b2 = EOS_Parameter.C2;
        amrex::Real KE = 0.5*( AMREX_D_TERM(u*u, +v*v, +w*w) );
        amrex::Real E1 = Conservative[INDEX_Energy1] - Alpha1_rho1*KE - Alpha1_rho1*D1;
        amrex::Real E2 = Conservative[INDEX_Energy2] - Alpha2_rho2*KE - Alpha2_rho2*D2;
        amrex::Real E  = E1 + E2;
        amrex::Real pinf1 = B1/(A1 + 1.0);
        amrex::Real pinf2 = B2/(A2 + 1.0);
        amrex::Real w1 = C1*Alpha1_rho1/A1;
        amrex::Real w2 = C2*Alpha2_rho2/A2;
        amrex::Real Phi = 1.0 - b1*Alpha1_rho1 - b2*Alpha2_rho2;

        amrex::Real aPT = Phi*(w1*A1 + w2*A2);
        amrex::Real bPT = Phi*(w1*(A1*pinf2 + B1) + w2*(A2*pinf1 + B2)) - E*(w1 + w2);
        amrex::Real cPT = Phi*(w1*B1*pinf2 + w2*B2*pinf1) - E*(w1*pinf2 + w2*pinf1);
        amrex::Real dPT = std::sqrt(bPT*bPT - 4.0*aPT*cPT);
        amrex::Real pR = (bPT > 0.0) ? -2.0*cPT/(bPT + dPT) : (-bPT + dPT)/(2.0*aPT);
        amrex::Real W1 = w1*(pR + pinf2);
        amrex::Real W2 = w2*(pR + pinf1);
        Alpha1R = b1*Alpha1_rho1 + Phi*W1/(W1 + W2);
        Alpha2R = b2*Alpha2_rho2 + Phi*W2/(W1 + W2);
        if (std::isnan(Alpha1R) && !std::isnan(Alpha2R)){
            Alpha1R = 1.0 - Alpha2R;
        }
        else if (!std::isnan(Alpha1R) && std::isnan(Alpha2R)){
            Alpha2R = 1.0 - Alpha1R;
        }
        amrex::Real Sum_AlphaR = Alpha1R + Alpha2R;
        Alpha1R /= Sum_AlphaR;
        Alpha2R /= Sum_AlphaR;
        amrex::Real Phi1R = NASG_Xi(Alpha1R, Alpha1_rho1, b1)*Alpha1R;
        amrex::Real Phi2R = NASG_Xi(Alpha2R, Alpha2_rho2, b2)*Alpha2R;
        pR = (E - B1*Phi1R - B2*Phi2R)/(A1*Phi1R + A2*Phi2R);
        P_Relaxed  = pR;
        Alpha1_p1R = Alpha1R*pR;
        Alpha2_p2R = Alpha2R*pR;
    }
    else{
        // Stiffened gas, ideal gas and Mie-Gruneisen (unchanged)
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
        amrex::Real dPT = std::sqrt(bPT*bPT - 4.0*aPT*cPT);
        amrex::Real pR = (bPT > 0.0) ? -2.0*cPT/(bPT + dPT) : (-bPT + dPT)/(2.0*aPT);
        Alpha1R = C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0))/(C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0)) + C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0)));
        Alpha2R = C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0))/(C1/A1*Alpha1_rho1*(pR + B2/(A2 + 1.0)) + C2/A2*Alpha2_rho2*(pR + B1/(A1 + 1.0)));
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
        P_Relaxed  = pR;
        Alpha1_p1R = Alpha1R*pR;
        Alpha2_p2R = Alpha2R*pR;
    }

    // The pressure from the energy above -pinf_k where alpha_k rho_k > 0, and a positive phase
    // temperature there, (alpha_k p_k + alpha_k pinf_k) > 0 (the two agree where alpha_kR > 0)
    amrex::Real pinf_1 = B1/(A1 + 1.0);
    amrex::Real pinf_2 = B2/(A2 + 1.0);
    bool admissible = (Alpha1_rho1 >= 0.0) && (Alpha2_rho2 >= 0.0)
                   && Relax_Finite(Alpha1R) && Relax_Finite(Alpha2R)
                   && Relax_Finite(Alpha1_p1R) && Relax_Finite(Alpha2_p2R);
    if (admissible && Alpha1_rho1 > 0.0) { admissible = (P_Relaxed + pinf_1 > 0.0) && (Alpha1_p1R + Alpha1R*pinf_1 > 0.0); }
    if (admissible && Alpha2_rho2 > 0.0) { admissible = (P_Relaxed + pinf_2 > 0.0) && (Alpha2_p2R + Alpha2R*pinf_2 > 0.0); }
    // NASG: below the covolume limit (not evaluated without a covolume)
    if (admissible && Alpha1_rho1 > 0.0 && EOS_Parameter.C1 > 0.0) { admissible = (Alpha1R > EOS_Parameter.C1*Alpha1_rho1); }
    if (admissible && Alpha2_rho2 > 0.0 && EOS_Parameter.C2 > 0.0) { admissible = (Alpha2R > EOS_Parameter.C2*Alpha2_rho2); }

    if (!admissible){
        Relaxation_Pressure_Instantaneous(ConservativeR, Conservative, EOS_Parameter, Iter, Path);
        Path += 10;
        return;
    }

    Iter = 0;
    Path = Relax_PathModel;
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
void Relaxation_PressureTemperature_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                                  amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                                  Parm const& EOS_Parameter)
{
    int Iter, Path;
    Relaxation_PressureTemperature_Instantaneous(ConservativeR, Conservative, EOS_Parameter, Iter, Path);
}

#endif
#endif