// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_ChemicalPotential_Adv_D_H_
#define FDM_ChemicalPotential_Adv_D_H_

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_ChemicalPotential_Adv_D(amrex::Real& ChemicalPotential,
                                 amrex::Real const& PF_c,
                                 amrex::Real const& PF_xm,
                                 amrex::Real const& PF_xp,
                                 amrex::Real const& Normalx_c,
                                 amrex::Real const& Normalx_xm,
                                 amrex::Real const& Normalx_xp,
#if (AMREX_SPACEDIM > 1)
                                 amrex::Real const& PF_ym,
                                 amrex::Real const& PF_yp,
                                 amrex::Real const& Normaly_c,
                                 amrex::Real const& Normaly_ym,
                                 amrex::Real const& Normaly_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                                 amrex::Real const& PF_zm,
                                 amrex::Real const& PF_zp,
                                 amrex::Real const& Normalz_c,
                                 amrex::Real const& Normalz_zm,
                                 amrex::Real const& Normalz_zp,
#endif
                                 amrex::Real const& Lambda,
                                 amrex::Real const& Eta,
                    AMREX_D_DECL(amrex::Real const& dx,
                                 amrex::Real const& dy,
                                 amrex::Real const& dz),
                                 int ID_Potential = 0)
{
    amrex::Real Laplace_PF;
    FDM_Laplace_D(Laplace_PF,
                  PF_c,
                  PF_xm, PF_xp,
#if (AMREX_SPACEDIM > 1)
                  PF_ym, PF_yp,
#endif
#if (AMREX_SPACEDIM > 2)
                  PF_zm, PF_zp,
#endif
     AMREX_D_DECL(dx, dy, dz));
    
    amrex::Real Normalx_fp;
    Interpolation2Face_Linear_D(Normalx_fp,
                                Normalx_c, Normalx_xp);
    amrex::Real Normalx_fm;
    Interpolation2Face_Linear_D(Normalx_fm,
                                Normalx_xm, Normalx_c);
#if (AMREX_SPACEDIM > 1)
    amrex::Real Normaly_fp;
    Interpolation2Face_Linear_D(Normaly_fp,
                                Normaly_c, Normaly_yp);
    amrex::Real Normaly_fm;
    Interpolation2Face_Linear_D(Normaly_fm,
                                Normaly_ym, Normaly_c);
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Real Normalz_fp;
    Interpolation2Face_Linear_D(Normalz_fp,
                                Normalz_c, Normalz_zp);
    amrex::Real Normalz_fm;
    Interpolation2Face_Linear_D(Normalz_fm,
                                Normalz_zm, Normalz_c);
#endif
    amrex::Real Div_Normal;
    FVM_Divergence_D(Div_Normal,
                     Normalx_fp, Normalx_fm,
#if (AMREX_SPACEDIM > 1)
                     Normaly_fp, Normaly_fm,
#endif
#if (AMREX_SPACEDIM > 2)
                     Normalz_fp, Normalz_fm,
#endif
        AMREX_D_DECL(dx, dy, dz));
    amrex::Real Curvature_PF = -Div_Normal;
    
    amrex::Real Abs_Grad_PF = 0.0;
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
        amrex::Real PF_m;
        amrex::Real PF_p;
        amrex::Real dX;
        //x
        if (idim == 0){
            PF_m = PF_xm;
            PF_p = PF_xp;
            dX = dx;
        }
#if (AMREX_SPACEDIM > 1)
        //y
        else if (idim == 1){
            PF_m = PF_ym;
            PF_p = PF_yp;
            dX = dy;
        }
#endif
#if (AMREX_SPACEDIM > 2)
        //z
        else if (idim == 2){
            PF_m = PF_zm;
            PF_p = PF_zp;
            dX = dz;
        }
#endif
        amrex::Real Grad_PF;
        FDM_Gradient_CenterFromCenter_D(Grad_PF,
                                        PF_m, PF_c, PF_p,
                                        dX);
        Abs_Grad_PF += Grad_PF*Grad_PF;
    }
    Abs_Grad_PF = std::sqrt(Abs_Grad_PF);
    
    ChemicalPotential_Adv(ChemicalPotential,
                          PF_c, Laplace_PF, Abs_Grad_PF, Curvature_PF,
                          Lambda, Eta,
                          ID_Potential);
    
}

#endif