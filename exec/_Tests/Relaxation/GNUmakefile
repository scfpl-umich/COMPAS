# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Six-equation relaxations: 40 hard states per model, relaxed once (conservation, passive
# phases keeping their volume fraction, sum of the volume fractions, pressure equality,
# iterations) and stepped, for each equation of state. One build per model.
#
#   make              main2d.gnu.MPI.ex        two-phase six-equation model (BUILD=6Eq, the default)
#   make BUILD=6Eq-N  main2d.gnu.MPI.6Eq-N.ex  N-phase six-equation model, NPHASE = 3
#
# Each build has its own executable and object directory (USERSuffix), so the builds can sit
# side by side here. variants.txt lists the runs: the label, BUILD and inputs file of each
# base run, and its variants. For example
#   make -j4 && ./main2d.gnu.MPI.ex prob/inputs   (Relaxation-6Eq)
#   make -j4 BUILD=6Eq-N && ./main2d.gnu.MPI.6Eq-N.ex prob/inputs.Relaxation-6Eq-N
# "make BUILD=6Eq-N cleanconfig" removes one build, "make realclean" every build.

BUILDS := 6Eq 6Eq-N
BUILD  ?= $(firstword $(BUILDS))

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_HIP    = FALSE
USE_SYCL   = FALSE

# The default build keeps the plain executable name, the others add their BUILD
ifneq ($(BUILD),$(firstword $(BUILDS)))
USERSuffix := .$(BUILD)
endif

Bpack   := ./prob/Make.package
Blocs   := . ./prob

include ../Make.exec

#############################################
# CODE SPECIFIC DEFINITIONS
#############################################

ifeq ($(BUILD),6Eq)
DEFINES += -DPHYSICS=SIXEQS
else ifeq ($(BUILD),6Eq-N)
DEFINES += -DPHYSICS=SIXEQS_IE_NPHASE
DEFINES += -DNPHASE=3
else
$(error BUILD=$(BUILD) is not one of the builds of this case: $(BUILDS))
endif

DEFINES += -DADVECTION=true
DEFINES += -DDIFFUSION=false
DEFINES += -DNONCONSERVATIVE=true
