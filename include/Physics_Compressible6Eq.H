// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Compressible6Eq_H_
#define Physics_Compressible6Eq_H_

#include<EquationOfState_NASG.H>
#include<RelaxationPressure_Safeguarded.H>
#include<Variables_Compressible6Eq.H>
#include<Physics_Parameter_Compressible6Eq.H>

#include<Vector2Component_Compressible6Eq.H>
#include<Component2Vector_Compressible6Eq.H>
#include<VectorPhase2Vector_Compressible6Eq.H>

#include<Primitive2State_Compressible6Eq.H>
#include<State2Primitive_Compressible6Eq.H>
#include<Conservative2State_Compressible6Eq.H>
#include<State2Conservative_Compressible6Eq.H>
#include<EquationOfState_Phase_Compressible6Eq.H>
#include<Conservative2Primitive_Compressible6Eq.H>
#include<Primitive2Conservative_Compressible6Eq.H>
#include<Conservative2ConservativePhase_Compressible6Eq.H>
#include<ConservativePhase2Conservative_Compressible6Eq.H>
#include<Conservative2Flux_Compressible6Eq.H>
#include<rho_c2_Phase_Compressible6Eq.H>
#include<rho_c2_Mixture_Compressible6Eq.H>
#include<SoundSpeed_Compressible6Eq.H>
#include<Riemann_LLF_Compressible6Eq.H>
#include<Riemann_HLL_Compressible6Eq.H>
#include<Riemann_HLLC_Compressible6Eq.H>
#include<RiemannSolver_Compressible6Eq.H>

#include<Conservative2Velocity_Compressible6Eq.H>
#include<Conservative2Temperature_Compressible6Eq.H>
#include<Conservative2PhaseTemperature_Compressible6Eq.H>
#include<Conservative2Pressure_Compressible6Eq.H>
#include<Conservative2EnergyBound_Compressible6Eq.H>

#include<ComputeDerivedVar_Compressible6Eq.H>

#include<FluxDiffusion_Compressible6Eq.H>

#include<Relaxation_Pressure_Instantaneous_Compressible6Eq.H>
#include<Relaxation_PressureTemperature_Instantaneous_Compressible6Eq.H>
#include<Relaxation_Instantaneous_Compressible6Eq.H>


#endif