// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RiemannSolver_Compressible5Eq_H_
#define RiemannSolver_Compressible5Eq_H_

#include<Macros.H>
#if (PHYSICS==FIVEEQS)

template <typename TF, typename TV, typename TU>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void RiemannSolver(TF& Flux_Riemann,
                   TV& Vel_Riemann,
                   amrex::Real& cmax,
                   TU const& ConservativeL,
                   TU const& ConservativeR,
                   Parm const& EOS_Parameter,
                   int ID_Dimension,
                   int ID_K,
                   int ID_Riemann)
{
    //LLF
    if (ID_Riemann == 0){
        Riemann_LLF(Flux_Riemann, Vel_Riemann, cmax,
                    ConservativeL, ConservativeR,
                    EOS_Parameter,
                    ID_Dimension, ID_K);
    }
    //HLL
    else if (ID_Riemann == 1){
        Riemann_HLL(Flux_Riemann, Vel_Riemann, cmax,
                    ConservativeL, ConservativeR,
                    EOS_Parameter,
                    ID_Dimension, ID_K);
    }
    //HLLC
    else if (ID_Riemann == 2){
        Riemann_HLLC(Flux_Riemann, Vel_Riemann, cmax,
                     ConservativeL, ConservativeR,
                     EOS_Parameter,
                     ID_Dimension, ID_K);
    }
}

#endif
#endif