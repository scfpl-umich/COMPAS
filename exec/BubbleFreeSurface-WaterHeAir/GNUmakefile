# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

# DIM        = 3
DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE
USE_HIP    = FALSE
USE_SYCL   = FALSE

Bpack   := ./prob/Make.package 
Blocs   := . ./prob


include ../Make.exec


#############################################
# CODE SPECIFIC DEFINTIONS
#############################################

# N-phase six-equation model, internal-energy formulation. The same inputs also run with
# the five-equation N-phase model (PHYSICS=FIVEEQS_NPHASE, which ignores Physics.pressure_relaxation)
DEFINES += -DPHYSICS=SIXEQS_IE_NPHASE
# DEFINES += -DPHYSICS=FIVEEQS_NPHASE

# Water, bubble gas and air
DEFINES += -DNPHASE=3

DEFINES += -DADVECTION=true
DEFINES += -DDIFFUSION=false
DEFINES += -DNONCONSERVATIVE=true


