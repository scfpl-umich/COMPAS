# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

# Noble-Abel stiffened-gas EOS (prob.EOS.EOS = 3), one build per model: the air-water and
# water shock tubes, a water drop at uniform pressure and velocity (prob.problem = 1), and
# b = q = 0 against the stiffened gas.
#
#   make              main2d.gnu.MPI.ex        five-equation model (BUILD=5Eq, the default)
#   make BUILD=5Eq-N  main2d.gnu.MPI.5Eq-N.ex  N-phase five-equation model, NPHASE = 2
#   make BUILD=6Eq    main2d.gnu.MPI.6Eq.ex    two-phase six-equation model
#   make BUILD=6Eq-N  main2d.gnu.MPI.6Eq-N.ex  N-phase six-equation model, NPHASE = 2
#
# Each build has its own executable and object directory (USERSuffix), so the builds can sit
# side by side here. variants.txt lists the runs: the label, BUILD and inputs file of each
# base run, and its variants. For example
#   make -j4 && ./main2d.gnu.MPI.ex prob/inputs   (NASG-5Eq)
#   make -j4 BUILD=5Eq-N && ./main2d.gnu.MPI.5Eq-N.ex prob/inputs.NASG-5Eq-N
# "make BUILD=5Eq-N cleanconfig" removes one build, "make realclean" every build.

BUILDS := 5Eq 5Eq-N 6Eq 6Eq-N
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
