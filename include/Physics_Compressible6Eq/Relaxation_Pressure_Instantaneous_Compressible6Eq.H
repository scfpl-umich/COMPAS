// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Relaxation_Pressure_Instantaneous_Compressible6Eq_H_
#define Relaxation_Pressure_Instantaneous_Compressible6Eq_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS)

// Instantaneous pressure relaxation (P_I = the relaxed pressure), in closed form: a quadratic for
// the relaxed pressure, then the volume fractions, and the pressure again from the conservation of
// E_1 + E_2. The closed form is kept wherever it does not fail: a finite result, not 0/0 (the
// spurious root of an absent phase), and the root of the quadratic well conditioned: farther than
// 1e-8 (|P| + max pinf_k) from the largest asymptote -pinf_m of the present phases. Closer, the root
// loses most of its digits (cancellation in -b + sqrt(b^2 - 4ac) when the other root is near -pinf
// of the other phase, e.g. a gas near vacuum or a gas trace in a liquid in tension, and P + pinf_m
// not representable next to a Mie-Gruneisen asymptote), with volume fraction errors up to 0.1 and
// more. A volume fraction out of [0, 1] is not a failure: no volume fraction is clipped, and the
// bounds come from the transport. Where the closed form fails, the safeguarded relaxation of
// RelaxationPressure_Safeguarded.H gives the volume fractions, and the phase energies are reset at
// the common pressure from E_1 + E_2, with the masses and the momentum unchanged and E_1 + E_2
// conserved to round-off (the closed form conserves it to the round-off of its pressure update,
// about 1e-13 relative where B_k alpha_k dominates the phase energy). A phase with
// 0 < alpha_k <= Physics.pressure_relaxation_trace (default 0) also takes the safeguarded path,
// where it keeps its volume fraction. Iter and Path report the iterations and the path
// (Relax_Path*).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Pressure_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
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

    amrex::Real A2 = EOS_Parameter.A2;
    amrex::Real B2 = EOS_Parameter.B2;

    amrex::Real Alpha1R;
    amrex::Real Alpha2R;
    amrex::Real Alpha1_p1R;
    amrex::Real Alpha2_p2R;
    // Partial-form phi_k and Pi_k = phi_k p_k of each phase, for the safeguarded relaxation
    // (phi_k = alpha_k and Pi_k = alpha_k p_k without a covolume)
    amrex::Real Phi1_In = Alpha1;
    amrex::Real Phi2_In = Alpha2;
    amrex::Real Pi1_In  = Alpha1_p1;
    amrex::Real Pi2_In  = Alpha2_p2;
    // True where the closed form gave 0/0 for a volume fraction: an absent phase whose asymptote
    // P = -pinf_k the quadratic picked as its root
    bool RawNaN = false;
    // The root of the quadratic
    amrex::Real pQ;
    if (EOS_Parameter.C1 > 0.0 || EOS_Parameter.C2 > 0.0){
        // NASG (b_1 or b_2 > 0). With phi_k = alpha_k - b_k alpha_k rho_k, Pi_k = phi_k p_k and
        // E_k = alpha_k rho_k e_k - alpha_k rho_k q_k = A_k Pi_k + B_k phi_k, the energy of each phase
        // relaxed with the final pressure p, phi_kR (A_k p + B_k) - E_k = -p (alpha_kR - alpha_k), gives
        //   alpha_kR = (A_k Pi_k + B_k phi_k + p alpha_k + b_k alpha_k rho_k (A_k p + B_k))/((A_k+1) p + B_k).
        // The condition alpha_1R + alpha_2R = 1 gives the quadratic aP p^2 + bP p + cP = 0 below, which is
        // the stiffened-gas one with alpha_k -> phi_k, alpha_k p_k -> Pi_k and A1 A2 -> A1 A2 (phi_1 + phi_2);
        // the same root is taken. The pressure is then recomputed from the conservation of E_1 + E_2.
        // E_k and phi_k come from the conservative state without the covolume floor, so that a phase
        // compressed past 1/b_k in a mixture cell by the hyperbolic step (phi_k < 0) is relaxed back.
        amrex::Real b1 = EOS_Parameter.C1;
        amrex::Real b2 = EOS_Parameter.C2;
        amrex::Real KE = 0.5*( AMREX_D_TERM(u*u, +v*v, +w*w) );
        amrex::Real E1 = Conservative[INDEX_Energy1] - Alpha1_rho1*KE - Alpha1_rho1*EOS_Parameter.D1;
        amrex::Real E2 = Conservative[INDEX_Energy2] - Alpha2_rho2*KE - Alpha2_rho2*EOS_Parameter.D2;
        amrex::Real Phi1 = Alpha1 - b1*Alpha1_rho1;
        amrex::Real Phi2 = Alpha2 - b2*Alpha2_rho2;
        amrex::Real Pi1  = (E1 - B1*Phi1)/A1;
        amrex::Real Pi2  = (E2 - B2*Phi2)/A2;
        Phi1_In = Phi1;
        Phi2_In = Phi2;
        Pi1_In  = Pi1;
        Pi2_In  = Pi2;

        amrex::Real aP = A1*A2*(Phi1 + Phi2) + A1*Phi1 + A2*Phi2;
        amrex::Real bP = A1*B2*Phi1 + A2*B1*Phi2 - (A2 + 1.0)*A1*Pi1 - (A1 + 1.0)*A2*Pi2;
        amrex::Real cP = -B2*A1*Pi1 - B1*A2*Pi2;
        amrex::Real pR = (-bP + std::sqrt(bP*bP - 4.0*aP*cP))/(2.0*aP);
        pQ = pR;
        Alpha1R = (E1 + pR*Alpha1 + b1*Alpha1_rho1*(A1*pR + B1))/((A1 + 1.0)*pR + B1);
        Alpha2R = (E2 + pR*Alpha2 + b2*Alpha2_rho2*(A2*pR + B2))/((A2 + 1.0)*pR + B2);
        RawNaN = std::isnan(Alpha1R) || std::isnan(Alpha2R);
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
        pR = (E1 + E2 - B1*Phi1R - B2*Phi2R)/(A1*Phi1R + A2*Phi2R);
        Alpha1_p1R = Alpha1R*pR;
        Alpha2_p2R = Alpha2R*pR;
    }
    else{
        // Stiffened gas, ideal gas and Mie-Gruneisen (unchanged)
        amrex::Real aP = A1*A2 + A1*Alpha1 + A2*Alpha2;
        amrex::Real bP = A1*B2*Alpha1 + A2*B1*Alpha2 - (A2 + 1.0)*A1*Alpha1_p1 - (A1 + 1.0)*A2*Alpha2_p2;
        amrex::Real cP = -B2*A1*Alpha1_p1 - B1*A2*Alpha2_p2;
        amrex::Real pR = (-bP + std::sqrt(bP*bP - 4.0*aP*cP))/(2.0*aP);
        pQ = pR;
        Alpha1R = (A1*Alpha1_p1 + B1*Alpha1 + pR*Alpha1)/((A1 + 1.0)*pR + B1);
        Alpha2R = (A2*Alpha2_p2 + B2*Alpha2 + pR*Alpha2)/((A2 + 1.0)*pR + B2);
        RawNaN = std::isnan(Alpha1R) || std::isnan(Alpha2R);
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
        Alpha1_p1R = Alpha1R*pR;
        Alpha2_p2R = Alpha2R*pR;
    }

    // With an absent phase the quadratic also has the root P = -pinf_k of that phase, which the
    // closed form picks when it lies above the true root (possible for a Mie-Gruneisen phase with
    // p_ref > 0). The normalization and the pressure from the energy repair that, unless the
    // volume fraction of the absent phase is exactly 0/0 (RawNaN): the old patch 1 - alpha_other
    // would then give it the volume
    amrex::Real Alpha_Trace = EOS_Parameter.Physics_Parm.pressure_relaxation_trace;
    // Conditioning of the root: distance to the largest asymptote of the present phases
    amrex::Real pinf_1 = B1/(A1 + 1.0);
    amrex::Real pinf_2 = B2/(A2 + 1.0);
    amrex::Real pinf_m = std::numeric_limits<amrex::Real>::max();
    amrex::Real pinf_scale = 0.0;
    if (Alpha1 > 0.0) { pinf_m = amrex::min(pinf_m, pinf_1); pinf_scale = amrex::max(pinf_scale, std::abs(pinf_1)); }
    if (Alpha2 > 0.0) { pinf_m = amrex::min(pinf_m, pinf_2); pinf_scale = amrex::max(pinf_scale, std::abs(pinf_2)); }
    bool well_conditioned = (std::abs(pQ + pinf_m) > 1.0e-8*(std::abs(pQ) + pinf_scale));
    bool admissible = !RawNaN && well_conditioned
                   && Relax_Finite(Alpha1R) && Relax_Finite(Alpha1_p1R) && Relax_Finite(Alpha2_p2R)
                   && !(Alpha1 > 0.0 && Alpha1 <= Alpha_Trace) && !(Alpha2 > 0.0 && Alpha2 <= Alpha_Trace);

    if (admissible){
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
        return;
    }

    // Safeguarded relaxation
    amrex::Array<amrex::Real,2> AlphaR;
    amrex::Real PR;
    Path = RelaxPressure_Safeguarded<2>(AlphaR, PR, Iter,
                                        {Alpha1, Alpha2}, {Alpha1_rho1, Alpha2_rho2},
                                        {Phi1_In, Phi2_In}, {Pi1_In, Pi2_In},
                                        {A1, A2}, {B1, B2}, Alpha_Trace);

    // Phase energies at the common pressure from E_1 + E_2; masses and momentum unchanged. The state
    // holds alpha_1 only, so phase 2 takes alpha_2 = 1 - alpha_1, as everywhere else. Unless both
    // phases are active (Relax_PathSolve) the volume fractions do not change: a passive phase keeps
    // its alpha_k, and the other one takes 1 - alpha_k, which is its own alpha in exact arithmetic;
    // alpha_1 is kept as stored, so that alpha_2 = 1 - alpha_1 is kept to the bit too
    if (Path != Relax_PathSolve) { AlphaR[0] = Alpha1; }
    AlphaR[1] = 1.0 - AlphaR[0];
    amrex::Real rho  = Alpha1_rho1 + Alpha2_rho2;
    amrex::Real rho_KE = 0.0;
    if (rho > 0.0){
        rho_KE = 0.5*( AMREX_D_TERM( Conservative[INDEX_MomentumX]*Conservative[INDEX_MomentumX],
                                    +Conservative[INDEX_MomentumY]*Conservative[INDEX_MomentumY],
                                    +Conservative[INDEX_MomentumZ]*Conservative[INDEX_MomentumZ] ) )/rho;
    }
    amrex::Real KE = (rho > 0.0) ? rho_KE/rho : 0.0;
    amrex::Real E  = Conservative[INDEX_Energy1] + Conservative[INDEX_Energy2];
    amrex::Real Phi1R = NASG_Xi(AlphaR[0], Alpha1_rho1, EOS_Parameter.C1)*AlphaR[0];
    amrex::Real Phi2R = NASG_Xi(AlphaR[1], Alpha2_rho2, EOS_Parameter.C2)*AlphaR[1];
    amrex::Real pR = (E - rho_KE - Alpha1_rho1*EOS_Parameter.D1 - Alpha2_rho2*EOS_Parameter.D2
                      - B1*Phi1R - B2*Phi2R)/(A1*Phi1R + A2*Phi2R);
    amrex::Real E1R = A1*Phi1R*pR + B1*Phi1R + Alpha1_rho1*(EOS_Parameter.D1 + KE);
    amrex::Real E2R = A2*Phi2R*pR + B2*Phi2R + Alpha2_rho2*(EOS_Parameter.D2 + KE);
    // The phase with the larger volume takes the rounding, so that E_1 + E_2 is kept
    if (AlphaR[0] >= AlphaR[1]) { E1R = E - E2R; }
    else                        { E2R = E - E1R; }

    ConservativeR = Conservative;
    ConservativeR[INDEX_Energy1]         = E1R;
    ConservativeR[INDEX_Energy2]         = E2R;
    ConservativeR[INDEX_VolumeFraction1] = AlphaR[0];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Relaxation_Pressure_Instantaneous(amrex::Array<amrex::Real,5+AMREX_SPACEDIM>& ConservativeR,
                                       amrex::Array<amrex::Real,5+AMREX_SPACEDIM> const& Conservative,
                                       Parm const& EOS_Parameter)
{
    int Iter, Path;
    Relaxation_Pressure_Instantaneous(ConservativeR, Conservative, EOS_Parameter, Iter, Path);
}

#endif
#endif