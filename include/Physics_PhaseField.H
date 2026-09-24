// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_PhaseField_H_
#define Physics_PhaseField_H_

#include<PhaseField_Parameter.H>

#include<Mapping_Forward.H>
#include<Mapping_Backward.H>

#include<PotentialFunction_DoubleWell.H>
#include<PotentialFunction_PhaseField.H>

#include<WeightFunction_Quadratic.H>
#include<WeightFunction_PhaseField.H>

#include<ChemicalPotential_VanDerWaals.H>
#include<ChemicalPotential_Adv.H>

#include<PhaseFieldMechanism_AC.H>
#include<PhaseFieldMechanism_CAC.H>

#include<PhaseFieldFlux_CDI.H>
#include<PhaseFieldFlux_ACDI.H>

#include<PhaseFieldFlux2GradQPF.H>

#endif