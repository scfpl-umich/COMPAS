// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_Conservative2FluxDiffusion_D_Compressible5Eq_NPhase_H_
#define FDM_Conservative2FluxDiffusion_D_Compressible5Eq_NPhase_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_Conservative2FluxDiffusion_D(amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1>& Flux,
                                      amrex::Real& dmax,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_c,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xm,
#if (AMREX_SPACEDIM > 1)
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_ym,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_yp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xm_ym,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xm_yp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xp_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_zm,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_zp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xm_zm,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xm_zp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_xp_zm,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_ym_zm,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_ym_zp,
                                      amrex::Array<amrex::Real,2*NPHASE+AMREX_SPACEDIM+1> const& Conservative_yp_zm,
#endif
                                      amrex::Real const& mu,
                                      amrex::Real const& mub,
                                      amrex::Real const& Kappa,
                         AMREX_D_DECL(amrex::Real const& dx,
                                      amrex::Real const& dy,
                                      amrex::Real const& dz),
                                      Parm const& EOS_Parameter,
                                      int ID_Dimension)
{
    //Conservative_c
    AMREX_D_TERM( amrex::Real u_c;,
                  amrex::Real v_c;,
                  amrex::Real w_c );
    Conservative2Velocity(AMREX_D_DECL(u_c, v_c, w_c),
                          Conservative_c);
    amrex::Real T_c;
    Conservative2Temperature(T_c,
                             Conservative_c,
                             EOS_Parameter);
    //Conservative_xm
    AMREX_D_TERM( amrex::Real u_xm;,
                  amrex::Real v_xm;,
                  amrex::Real w_xm );
    Conservative2Velocity(AMREX_D_DECL(u_xm, v_xm, w_xm),
                          Conservative_xm);
    amrex::Real T_xm;
    Conservative2Temperature(T_xm,
                             Conservative_xm,
                             EOS_Parameter);
#if (AMREX_SPACEDIM > 1)
    //Conservative_xp
    AMREX_D_TERM( amrex::Real u_xp;,
                  amrex::Real v_xp;,
                  amrex::Real w_xp );
    Conservative2Velocity(AMREX_D_DECL(u_xp, v_xp, w_xp),
                          Conservative_xp);
    amrex::Real T_xp;
    Conservative2Temperature(T_xp,
                             Conservative_xp,
                             EOS_Parameter);
    //Conservative_ym
    AMREX_D_TERM( amrex::Real u_ym;,
                  amrex::Real v_ym;,
                  amrex::Real w_ym );
    Conservative2Velocity(AMREX_D_DECL(u_ym, v_ym, w_ym),
                          Conservative_ym);
    amrex::Real T_ym;
    Conservative2Temperature(T_ym,
                             Conservative_ym,
                             EOS_Parameter);
    //Conservative_yp
    AMREX_D_TERM( amrex::Real u_yp;,
                  amrex::Real v_yp;,
                  amrex::Real w_yp );
    Conservative2Velocity(AMREX_D_DECL(u_yp, v_yp, w_yp),
                          Conservative_yp);
    amrex::Real T_yp;
    Conservative2Temperature(T_yp,
                             Conservative_yp,
                             EOS_Parameter);
    //Conservative_xm_ym
    AMREX_D_TERM( amrex::Real u_xm_ym;,
                  amrex::Real v_xm_ym;,
                  amrex::Real w_xm_ym );
    Conservative2Velocity(AMREX_D_DECL(u_xm_ym, v_xm_ym, w_xm_ym),
                          Conservative_xm_ym);
    amrex::Real T_xm_ym;
    Conservative2Temperature(T_xm_ym,
                             Conservative_xm_ym,
                             EOS_Parameter);
    //Conservative_xm_yp
    AMREX_D_TERM( amrex::Real u_xm_yp;,
                  amrex::Real v_xm_yp;,
                  amrex::Real w_xm_yp );
    Conservative2Velocity(AMREX_D_DECL(u_xm_yp, v_xm_yp, w_xm_yp),
                          Conservative_xm_yp);
    amrex::Real T_xm_yp;
    Conservative2Temperature(T_xm_yp,
                             Conservative_xm_yp,
                             EOS_Parameter);
    //Conservative_xp_ym
    AMREX_D_TERM( amrex::Real u_xp_ym;,
                  amrex::Real v_xp_ym;,
                  amrex::Real w_xp_ym );
    Conservative2Velocity(AMREX_D_DECL(u_xp_ym, v_xp_ym, w_xp_ym),
                          Conservative_xp_ym);
    amrex::Real T_xp_ym;
    Conservative2Temperature(T_xp_ym,
                             Conservative_xp_ym,
                             EOS_Parameter);
#endif
#if (AMREX_SPACEDIM > 2)
    //Conservative_zm
    AMREX_D_TERM( amrex::Real u_zm;,
                  amrex::Real v_zm;,
                  amrex::Real w_zm );
    Conservative2Velocity(AMREX_D_DECL(u_zm, v_zm, w_zm),
                          Conservative_zm);
    amrex::Real T_zm;
    Conservative2Temperature(T_zm,
                             Conservative_zm,
                             EOS_Parameter);
    //Conservative_zp
    AMREX_D_TERM( amrex::Real u_zp;,
                  amrex::Real v_zp;,
                  amrex::Real w_zp );
    Conservative2Velocity(AMREX_D_DECL(u_zp, v_zp, w_zp),
                          Conservative_zp);
    amrex::Real T_zp;
    Conservative2Temperature(T_zp,
                             Conservative_zp,
                             EOS_Parameter);
    //Conservative_xm_zm
    AMREX_D_TERM( amrex::Real u_xm_zm;,
                  amrex::Real v_xm_zm;,
                  amrex::Real w_xm_zm );
    Conservative2Velocity(AMREX_D_DECL(u_xm_zm, v_xm_zm, w_xm_zm),
                          Conservative_xm_zm);
    amrex::Real T_xm_zm;
    Conservative2Temperature(T_xm_zm,
                             Conservative_xm_zm,
                             EOS_Parameter);
    //Conservative_xm_zp
    AMREX_D_TERM( amrex::Real u_xm_zp;,
                  amrex::Real v_xm_zp;,
                  amrex::Real w_xm_zp );
    Conservative2Velocity(AMREX_D_DECL(u_xm_zp, v_xm_zp, w_xm_zp),
                          Conservative_xm_zp);
    amrex::Real T_xm_zp;
    Conservative2Temperature(T_xm_zp,
                             Conservative_xm_zp,
                             EOS_Parameter);
    //Conservative_xp_zm
    AMREX_D_TERM( amrex::Real u_xp_zm;,
                  amrex::Real v_xp_zm;,
                  amrex::Real w_xp_zm );
    Conservative2Velocity(AMREX_D_DECL(u_xp_zm, v_xp_zm, w_xp_zm),
                          Conservative_xp_zm);
    amrex::Real T_xp_zm;
    Conservative2Temperature(T_xp_zm,
                             Conservative_xp_zm,
                             EOS_Parameter);
    //Conservative_ym_zm
    AMREX_D_TERM( amrex::Real u_ym_zm;,
                  amrex::Real v_ym_zm;,
                  amrex::Real w_ym_zm );
    Conservative2Velocity(AMREX_D_DECL(u_ym_zm, v_ym_zm, w_ym_zm),
                          Conservative_ym_zm);
    amrex::Real T_ym_zm;
    Conservative2Temperature(T_ym_zm,
                             Conservative_ym_zm,
                             EOS_Parameter);
    //Conservative_ym_zp
    AMREX_D_TERM( amrex::Real u_ym_zp;,
                  amrex::Real v_ym_zp;,
                  amrex::Real w_ym_zp );
    Conservative2Velocity(AMREX_D_DECL(u_ym_zp, v_ym_zp, w_ym_zp),
                          Conservative_ym_zp);
    amrex::Real T_ym_zp;
    Conservative2Temperature(T_ym_zp,
                             Conservative_ym_zp,
                             EOS_Parameter);
    //Conservative_yp_zm
    AMREX_D_TERM( amrex::Real u_yp_zm;,
                  amrex::Real v_yp_zm;,
                  amrex::Real w_yp_zm );
    Conservative2Velocity(AMREX_D_DECL(u_yp_zm, v_yp_zm, w_yp_zm),
                          Conservative_yp_zm);
    amrex::Real T_yp_zm;
    Conservative2Temperature(T_yp_zm,
                             Conservative_yp_zm,
                             EOS_Parameter);
#endif
    
    FDM_FluxDiffusion_D(Flux,
                        dmax,
                        u_c,
                        u_xm,
                        T_c,
                        T_xm,
#if (AMREX_SPACEDIM > 1)
                        u_xp, u_ym, u_yp,
                        u_xm_ym, u_xm_yp, u_xp_ym,
                        v_c,
                        v_xm, v_xp, v_ym, v_yp,
                        v_xm_ym, v_xm_yp, v_xp_ym,
                        T_ym,
#endif
#if (AMREX_SPACEDIM > 2)
                        u_zm, u_zp,
                        u_xm_zm, u_xm_zp, u_xp_zm,
                        v_zm, v_zp,
                        v_ym_zm, v_ym_zp, v_yp_zm,
                        w_c,
                        w_xm, w_xp, w_ym, w_yp, w_zm, w_zp,
                        w_xm_zm, w_xm_zp, w_xp_zm,
                        w_ym_zm, w_ym_zp, w_yp_zm,
                        T_zm,
#endif
                        mu, mub, Kappa,
           AMREX_D_DECL(dx, dy, dz),
                        ID_Dimension);
    
}

#endif
#endif