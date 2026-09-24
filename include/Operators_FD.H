// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Operators_FD_H_
#define Operators_FD_H_

#include<Interpolation2Face_Linear_D.H>
#include<Interpolation2Face_Quick_D.H>
#include<Interpolation2Face_WENO3_D.H>
#include<Interpolation2Face_WENO5_D.H>
#include<Interpolation2Face_WENO3B_D.H>
#include<Interpolation2Face_WENO5B_D.H>

#include<Interpolation2Face_Linear_K.H>
#include<Interpolation2Face_Quick_K.H>
#include<Interpolation2Face_WENO3_K.H>
#include<Interpolation2Face_WENO5_K.H>
#include<Interpolation2Face_WENO3B_K.H>
#include<Interpolation2Face_WENO5B_K.H>

#include<Interpolation2Face_D.H>
#include<Interpolation2Face_K.H>

#include<Interpolation2Edge_Linear_D.H>
#include<Interpolation2Vertex_Linear_D.H>

#include<FDM_Gradient_D.H>

#include<FDM_GradientX_XFaceFromCenter_D.H>
#include<FDM_GradientY_XFaceFromCenter_D.H>
#include<FDM_GradientZ_XFaceFromCenter_D.H>
#include<FDM_Gradient_XFaceFromCenter_D.H>

#include<FDM_GradientX_YFaceFromCenter_D.H>
#include<FDM_GradientY_YFaceFromCenter_D.H>
#include<FDM_GradientZ_YFaceFromCenter_D.H>
#include<FDM_Gradient_YFaceFromCenter_D.H>

#include<FDM_GradientX_ZFaceFromCenter_D.H>
#include<FDM_GradientY_ZFaceFromCenter_D.H>
#include<FDM_GradientZ_ZFaceFromCenter_D.H>
#include<FDM_Gradient_ZFaceFromCenter_D.H>

#include<FDM_Gradient_CenterFromCenter_D.H>
#include<FDM_Gradient_CenterFromCenter_K.H>

#include<FDM_Gradient_FaceFromCenter_K.H>
#include<FDM_Gradient_XFaceFromCenter_K.H>
#include<FDM_Gradient_YFaceFromCenter_K.H>
#include<FDM_Gradient_ZFaceFromCenter_K.H>

#include<FDM_Normal_CenterFromCenter_D.H>
#include<FDM_Normal_XFaceFromCenter_D.H>
#include<FDM_Normal_YFaceFromCenter_D.H>
#include<FDM_Normal_ZFaceFromCenter_D.H>
#include<FDM_Curvature_CenterFromCenter_D.H>

#include<FDM_Normal_CenterFromCenter_K.H>
#include<FDM_Normal_XFaceFromCenter_K.H>
#include<FDM_Normal_YFaceFromCenter_K.H>
#include<FDM_Normal_ZFaceFromCenter_K.H>
#include<FDM_Curvature_CenterFromCenter_K.H>

#include<FDM_Laplace_D.H>
#include<FDM_Laplace_K.H>

#include<THINC_Select_K.H>

#endif