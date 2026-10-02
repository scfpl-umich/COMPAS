// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Compressible5Eq_H_
#define Physics_Compressible5Eq_H_

#include<EquationOfState_NASG.H>
#include<Variables_Compressible5Eq.H>
#include<Physics_Parameter_Compressible5Eq.H>

#include<Vector2Component_Compressible5Eq.H>
#include<Component2Vector_Compressible5Eq.H>
#include<VectorPhase2Vector_Compressible5Eq.H>

#include<Primitive2State_Compressible5Eq.H>
#include<State2Primitive_Compressible5Eq.H>
#include<Conservative2State_Compressible5Eq.H>
#include<State2Conservative_Compressible5Eq.H>
#include<EquationOfState_Mixture_Compressible5Eq.H>
#include<Conservative2Primitive_Compressible5Eq.H>
#include<Primitive2Conservative_Compressible5Eq.H>
#include<Conservative2Flux_Compressible5Eq.H>
#include<rho_c2_Phase_Compressible5Eq.H>
#include<rho_c2_Mixture_Compressible5Eq.H>
#include<SoundSpeed_Compressible5Eq.H>
#include<NonConservative_Compressible5Eq.H>
#include<Riemann_LLF_Compressible5Eq.H>
#include<Riemann_HLL_Compressible5Eq.H>
#include<Riemann_HLLC_Compressible5Eq.H>
#include<Primitive2CharacteristicDecomposition_Compressible5Eq.H>
#include<AnalyticHyperbolicFlux_Compressible5Eq.H>
#include<RiemannSolver_Compressible5Eq.H>
#include<ComputeDerivedVar_Compressible5Eq.H>

#include<Conservative2Velocity_Compressible5Eq.H>
#include<Conservative2Temperature_Compressible5Eq.H>
#include<Conservative2Pressure_Compressible5Eq.H>
#include<FluxDiffusion_Compressible5Eq.H>
#include<ConservativeGradients2PrimitiveGradients_Compressible5Eq_NPhase.H>
#include<AnalyticDiffusiveFlux_Compressible5Eq.H>

#include<EquationOfState_Phase_Compressible5Eq.H>
#include<Conservative2ConservativePhase_Compressible5Eq.H>
#include<ConservativePhase2Conservative_Compressible5Eq.H>

#include<Conservative2A_Compressible5Eq.H>
#include<Conservative2Pmin_Compressible5Eq.H>
#include<Conservative2EnergyBound_Compressible5Eq.H>
#include<Conservative2Lambda_HLL_Compressible5Eq.H>
#include<Riemann_HLL_Bound_Compressible5Eq.H>
#include<Conservative2K_Compressible5Eq.H>

#include<RefluxSourceTerms_Compressible5Eq.H>
#include<AnalyticNonConservativeFlux_Compressible5Eq.H>

#endif