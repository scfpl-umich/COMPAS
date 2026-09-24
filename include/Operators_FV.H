// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Operators_FV_H_
#define Operators_FV_H_

#include<FiniteVolume_Parameter.H>

#include<Quadrature.H>

#include<f_Downwind.H>
#include<f_Upwind.H>
#include<f_LimiterBound.H>
#include<fC_QuadratureRule.H>

#include<Limiter_Minmod.H>
#include<Limiter_vanLeer.H>
#include<Limiter_Superbee.H>
#include<Limiter_MC.H>
#include<Slope_TVD.H>

#include<Limiter_BoundMax.H>
#include<Limiter_BoundMin.H>
#include<Limiter_BoundMaxMin.H>
#include<Limiter_Bound.H>
#include<Limiter_Consistent.H>

#include<WENO3.H>
#include<WENO5.H>

#include<Reconstruction2Face_Godunov.H>
#include<Reconstruction2Face_Linear.H>
#include<Reconstruction2Face_WENO.H>

#include<Reconstruction_Godunov.H>
#include<Reconstruction_Linear.H>
#include<Reconstruction_CenterFace2.H>
#include<Reconstruction_CenterFace4.H>

#include<Reconstruction2Face_Godunov_D.H>
#include<Reconstruction2Face_MUSCL_D.H>
#include<Reconstruction2Face_WENO3_D.H>
#include<Reconstruction2Face_WENO5_D.H>
#include<Reconstruction2Face_WENO3B_D.H>
#include<Reconstruction2Face_WENO5B_D.H>

#include<Reconstruction2Face_Godunov_K.H>
#include<Reconstruction2Face_MUSCL_K.H>
#include<Reconstruction2Face_WENO3_K.H>
#include<Reconstruction2Face_WENO5_K.H>
#include<Reconstruction2Face_WENO3B_K.H>
#include<Reconstruction2Face_WENO5B_K.H>
#include<Reconstruction2Face_THINC_D.H>
#include<Reconstruction2Face_MTHINC_D.H>

#include<Reconstruction_Godunov_D.H>
#include<Reconstruction_MUSCL_D.H>
#include<Reconstruction_WENO3_D.H>
#include<Reconstruction_WENO5_D.H>
#include<Reconstruction_WENO3B_D.H>
#include<Reconstruction_WENO5B_D.H>

#include<Reconstruction_Godunov_K.H>
#include<Reconstruction_MUSCL_K.H>
#include<Reconstruction_WENO3_K.H>
#include<Reconstruction_WENO5_K.H>
#include<Reconstruction_WENO3B_K.H>
#include<Reconstruction_WENO5B_K.H>

#include<Reconstruction2Face_D.H>
#include<Reconstruction2Face_K.H>
#include<Reconstruction_D.H>
#include<Reconstruction_K.H>

#include<FVM_Divergence_D.H>
#include<FVM_Divergence_K.H>

#include<GetDifferentialMeasures.H>
#include<FVM_SurfaceIntegral.H>

#endif