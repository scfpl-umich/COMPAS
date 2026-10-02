// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Compressible6Eq_IE_NPhase_H_
#define Physics_Compressible6Eq_IE_NPhase_H_

#include<EquationOfState_NASG.H>
#include<Variables_Compressible6Eq_IE_NPhase.H>
#include<Physics_Parameter_Compressible6Eq_IE_NPhase.H>

#include<Vector2Component_Compressible6Eq_IE_NPhase.H>
#include<Component2Vector_Compressible6Eq_IE_NPhase.H>
#include<VectorPhase2Vector_Compressible6Eq_IE_NPhase.H>

#include<Primitive2State_Compressible6Eq_IE_NPhase.H>
#include<State2Primitive_Compressible6Eq_IE_NPhase.H>
#include<Conservative2State_Compressible6Eq_IE_NPhase.H>
#include<State2Conservative_Compressible6Eq_IE_NPhase.H>
#include<EquationOfState_Phase_Compressible6Eq_IE_NPhase.H>
#include<Conservative2Primitive_Compressible6Eq_IE_NPhase.H>
#include<Primitive2Conservative_Compressible6Eq_IE_NPhase.H>
#include<Conservative2ConservativePhase_Compressible6Eq_IE_NPhase.H>
#include<ConservativePhase2Conservative_Compressible6Eq_IE_NPhase.H>
#include<Conservative2Flux_Compressible6Eq_IE_NPhase.H>
#include<rho_c2_Phase_Compressible6Eq_IE_NPhase.H>
#include<rho_c2_Mixture_Compressible6Eq_IE_NPhase.H>
#include<SoundSpeed_Compressible6Eq_IE_NPhase.H>
#include<Riemann_LLF_Compressible6Eq_IE_NPhase.H>
#include<Riemann_HLL_Compressible6Eq_IE_NPhase.H>
#include<Riemann_HLLC_Compressible6Eq_IE_NPhase.H>
#include<RiemannSolver_Compressible6Eq_IE_NPhase.H>

#include<Conservative2Velocity_Compressible6Eq_IE_NPhase.H>
#include<Conservative2Temperature_Compressible6Eq_IE_NPhase.H>
#include<Conservative2PhaseTemperature_Compressible6Eq_IE_NPhase.H>
#include<Conservative2Pressure_Compressible6Eq_IE_NPhase.H>
#include<Conservative2EnergyBound_Compressible6Eq_IE_NPhase.H>
#include<FluxDiffusion_Compressible6Eq_IE_NPhase.H>

#include<ComputeDerivedVar_Compressible6Eq_IE_NPhase.H>

#include<EnergyCorrection_Compressible6Eq_IE_NPhase.H>

#include<RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase.H>

#include<RefluxSourceTerms_Compressible6Eq_IE_NPhase.H>

#endif