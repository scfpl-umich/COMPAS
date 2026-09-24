// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_LimiterBound_Compressible5Eq_H_
#define FVM_Reconstruction2Face_Conservative_LimiterBound_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_LimiterBound(amrex::Array<amrex::Real,4+AMREX_SPACEDIM>& Conservative_Face,
                                                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& Conservative_c,
                                                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeL_fp,
                                                       amrex::Array<amrex::Real,4+AMREX_SPACEDIM> const& ConservativeR_fm,
                                                       Parm const& EOS_Parameter,
                                                       amrex::Real const& Weight,
                                                       int ID_Face)
{
    amrex::Array<amrex::Real,4+AMREX_SPACEDIM> ConservativeC;
    for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
        amrex::Real f_c   = Conservative_c  [i_State];
        amrex::Real fL_fp = ConservativeL_fp[i_State];
        amrex::Real fR_fm = ConservativeR_fm[i_State];
        amrex::Real fC;
        fC_QuadratureRule(fC,
                          f_c, fL_fp, fR_fm,
                          Weight);
        ConservativeC[i_State] = fC;
    }

    amrex::Real EnergyBound_c;
    Conservative2EnergyBound(EnergyBound_c,
                             Conservative_c,
                             EOS_Parameter);
    amrex::Real EnergyBoundC;
    Conservative2EnergyBound(EnergyBoundC,
                             ConservativeC,
                             EOS_Parameter);
    amrex::Real EnergyBoundL_fp;
    Conservative2EnergyBound(EnergyBoundL_fp,
                             ConservativeL_fp,
                             EOS_Parameter);
    amrex::Real EnergyBoundR_fm;
    Conservative2EnergyBound(EnergyBoundR_fm,
                             ConservativeR_fm,
                             EOS_Parameter);
    amrex::Real EnergyBound_min = 0.0;
    amrex::Real ThetaB;
    Limiter_BoundMin(ThetaB,
                     EnergyBound_min,
                     EnergyBound_c, EnergyBoundC, EnergyBoundL_fp, EnergyBoundR_fm);
    
    //fL_fp
    if (ID_Face == 1){
        for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
            amrex::Real f   = ConservativeL_fp[i_State];
            amrex::Real f_c = Conservative_c  [i_State];
            amrex::Real f_Limit;
            f_LimiterBound(f_Limit,
                           f, f_c,
                           ThetaB);
            Conservative_Face[i_State] = f_Limit;
        }
    }
    //fR_fm
    else if (ID_Face == -1){
        for (int i_State = 0; i_State < 4+AMREX_SPACEDIM; i_State++){
            amrex::Real f   = ConservativeR_fm[i_State];
            amrex::Real f_c = Conservative_c  [i_State];
            amrex::Real f_Limit;
            f_LimiterBound(f_Limit,
                           f, f_c,
                           ThetaB);
            Conservative_Face[i_State] = f_Limit;
        }
    }
    
}

#endif
#endif