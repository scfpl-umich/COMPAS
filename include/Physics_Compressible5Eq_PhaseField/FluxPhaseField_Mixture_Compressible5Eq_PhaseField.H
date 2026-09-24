// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FluxPhaseField_Mixture_Compressible5Eq_PhaseField_H_
#define FluxPhaseField_Mixture_Compressible5Eq_PhaseField_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FluxPhaseField_Mixture(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Flux,
                            amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeL,
                            amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeR,
                            amrex::Real const& J1_rho1,
                            amrex::Real const& J2_rho2,
                            amrex::Real const& J1,
                            amrex::Real const& J2,
                            Parm const& EOS_Parameter)
{
    //Primitive
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> PrimitiveL;
    Conservative2Primitive(PrimitiveL,
                           ConservativeL,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1L;
    amrex::Real Alpha2_rho2L;
    AMREX_D_TERM( amrex::Real uL;,
                  amrex::Real vL;,
                  amrex::Real wL );
    amrex::Real pL;
    amrex::Real Alpha1L;
    amrex::Real Alpha2L;
    Primitive2State(Alpha1_rho1L,
                    Alpha2_rho2L,
       AMREX_D_DECL(uL, vL, wL),
                    pL,
                    Alpha1L,
                    PrimitiveL);
    Alpha2L = 1.0 - Alpha1L;
    
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> PrimitiveR;
    Conservative2Primitive(PrimitiveR,
                           ConservativeR,
                           EOS_Parameter);
    amrex::Real Alpha1_rho1R;
    amrex::Real Alpha2_rho2R;
    AMREX_D_TERM( amrex::Real uR;,
                  amrex::Real vR;,
                  amrex::Real wR );
    amrex::Real pR;
    amrex::Real Alpha1R;
    amrex::Real Alpha2R;
    Primitive2State(Alpha1_rho1R,
                    Alpha2_rho2R,
       AMREX_D_DECL(uR, vR, wR),
                    pR,
                    Alpha1R,
                    PrimitiveR);
    Alpha2R = 1.0 - Alpha1R;
    
    //Momentum
    amrex::Real M = J1_rho1 + J2_rho2;
    
    amrex::Real u_DM;
    f_Downwind(u_DM,
               uL,
               uR,
               M);
    amrex::Real M_u = M*u_DM;
#if (AMREX_SPACEDIM > 1)
    amrex::Real v_DM;
    f_Downwind(v_DM,
               vL,
               vR,
               M);
    amrex::Real M_v = M*v_DM;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real w_DM;
    f_Downwind(w_DM,
               wL,
               wR,
               M);
    amrex::Real M_w = M*w_DM;
#endif
    //Energy
    amrex::Real eK_DM = 0.5*( AMREX_D_TERM(u_DM*u_DM,+v_DM*v_DM,+w_DM*w_DM) );
    amrex::Real M_eK = M*eK_DM;
    
    amrex::Real p_DJ1;
    f_Downwind(p_DJ1,
               pL,
               pR,
               J1);
    amrex::Real J1_rho1_e1 = EOS_Parameter.A1*J1*p_DJ1 + EOS_Parameter.B1*J1 + EOS_Parameter.D1*J1_rho1;
    
    amrex::Real p_DJ2;
    f_Downwind(p_DJ2,
               pL,
               pR,
               J2);
    amrex::Real J2_rho2_e2 = EOS_Parameter.A2*J2*p_DJ2 + EOS_Parameter.B2*J2 + EOS_Parameter.D2*J2_rho2;
    
    //Flux_PhaseField
    Component2Vector(Flux,
                     J1_rho1,
                     J2_rho2,
        AMREX_D_DECL(M_u, M_v, M_w),
                     M_eK + J1_rho1_e1 + J2_rho2_e2,
                     J1);
    
}

#endif