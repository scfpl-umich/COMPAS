// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_SurfaceTension_H_
#define Physics_SurfaceTension_H_

// Surface tension kernels (SURFACE_TENSION = true), for every model. Include after the model's
// physics and operator headers.
#include<SurfaceTension_Parameter.H>
#include<FDM_Conservative2FluxCapillary_K.H>
#include<SurfaceTension_Pressure_K.H>

#endif
