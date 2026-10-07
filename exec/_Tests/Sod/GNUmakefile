# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Sod shock tube in 3D with AMR and outflow boundaries, one build per model; the variants
# run it along y and z, with walls, and the two-material water-air shock tube.
#
#   make              main3d.gnu.MPI.ex        five-equation model (BUILD=5Eq, the default)
#   make BUILD=5Eq-N  main3d.gnu.MPI.5Eq-N.ex  N-phase five-equation model, NPHASE = 2
#   make BUILD=6Eq    main3d.gnu.MPI.6Eq.ex    two-phase six-equation model
#   make BUILD=6Eq-N  main3d.gnu.MPI.6Eq-N.ex  N-phase six-equation model, NPHASE = 2
#
# Each build has its own executable and object directory (USERSuffix), so the builds can sit
# side by side here. variants.txt lists the runs: the label, BUILD and inputs file of each
# base run, and its variants. For example
#   make -j4 && ./main3d.gnu.MPI.ex prob/inputs   (Sod-5Eq)
#   make -j4 BUILD=5Eq-N && ./main3d.gnu.MPI.5Eq-N.ex prob/inputs.Sod-5Eq-N
# "make BUILD=5Eq-N cleanconfig" removes one build, "make realclean" every build.

BUILDS := 5Eq 5Eq-N 6Eq 6Eq-N
BUILD  ?= $(firstword $(BUILDS))

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 3

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

ifeq ($(BUILD),5Eq)
DEFINES += -DPHYSICS=FIVEEQS
else ifeq ($(BUILD),5Eq-N)
DEFINES += -DPHYSICS=FIVEEQS_NPHASE
DEFINES += -DNPHASE=2
else ifeq ($(BUILD),6Eq)
DEFINES += -DPHYSICS=SIXEQS
else ifeq ($(BUILD),6Eq-N)
DEFINES += -DPHYSICS=SIXEQS_IE_NPHASE
DEFINES += -DNPHASE=2
else
$(error BUILD=$(BUILD) is not one of the builds of this case: $(BUILDS))
endif

DEFINES += -DADVECTION=true
DEFINES += -DDIFFUSION=false
DEFINES += -DNONCONSERVATIVE=true
