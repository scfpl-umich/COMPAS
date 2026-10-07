// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Time step of all levels without subcycling (run.do_subcycle = 0), for the hyperbolic step and the
// Phase-Field step. The levels take every Runge-Kutta stage together: at each stage every level
// fills its ghost cells from the stage state (the coarse-fine ghost cells by spatial interpolation
// of the coarser level's stage state), computes its face fluxes and non-conservative face
// quantities, and with run.do_reflux = 1 the coarse faces covered by the finer level take the finer
// level's values averaged with the face areas, before any level is updated. The composite update is
// then conservative at every stage, and the non-conservative terms of a coarse cell use the same
// faces as its conservative fluxes; no reflux follows the step. Each stage ends with PostTimeStage
// (the pressure relaxation of the six-equation models) on every level and AverageDown. With one
// level the operations are those of the per-level integrators in the same order. The subcycled step
// (timeStepWithSubcycling) keeps the per-level integrators and the flux registers.

#include <COMPAS.H>
#include <Kernels.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Gpu.H>
#include <LinearAlgebra.H>
#include <Tools.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <Class_Compressible_PhaseField.H>

#include <Macros.H>

using namespace amrex;


void
Compressible_PhaseField::AdvanceAllLevels (Real time, Real dt_lev, int ncycle)
{
    AdvanceAllLevelsStageCoupled(dof_new, dof_old, nullptr, time, dt_lev, ncycle);
}


//Advance state vector at all levels with different RK schemes
void
Compressible_PhaseField::PhaseField_AdvanceAllLevels (Vector<MultiFab>& mf_new,
                                                      Vector<MultiFab>& mf_old,
                                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                                      Real time,
                                                      Real dt_lev,
                                                      int ncycle)
{
    AdvanceAllLevelsStageCoupled(mf_new, mf_old, &Grad_Q, time, dt_lev, ncycle);
}


// Energy-bound check of a stage state (FiniteVolume.ID_Bound = 1), as in the per-level integrators:
// true (and ID_Return = 1) if the step has to be restarted
bool
Compressible_PhaseField::StageBoundViolated (MultiFab const& U)
{
    if (h_parm->FiniteVolume_Parm.ID_Bound != 1){
        return false;
    }
    Real min_EnergyBound;
    MinEnergyBound(min_EnergyBound,
                   U,
                   *h_parm);
    h_parm->FiniteVolume_Parm.ID_Return = 0;
    if (min_EnergyBound < 0){
        h_parm->FiniteVolume_Parm.ID_Return = 1;
        amrex::Real min_Pressure;
        MinPressure(min_Pressure,
                    U,
                    *h_parm);
        amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
        amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
        return true;
    }
    return false;
}


// Right-hand side dUdt[lev] of every level for the stage state U_stage at stage_time (Uborder[lev]
// receives the stage state with its ghost cells). Grad_Q == nullptr: the hyperbolic step;
// otherwise the Phase-Field step with these face gradients. dt_stage (the stage weight times dt) is
// what the per-level integrators pass to the user source term and the CFL check
void
Compressible_PhaseField::StageRHSAllLevels (Vector<MultiFab>& dUdt,
                                            Vector<MultiFab>& Uborder,
                                            Vector<MultiFab>& U_stage,
                                            Vector<Array<MultiFab,AMREX_SPACEDIM>> const* Grad_Q,
                                            Real stage_time,
                                            Real dt_stage,
                                            int ncycle,
                                            int stage)
{
    BL_PROFILE("StageRHSAllLevels()");

    const int nlev = finest_level + 1;

    Vector<Array<MultiFab,AMREX_SPACEDIM>> fluxes(nlev);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> fluxes_nc(nlev);

    // 1. ghost cells from the stage state, 2. face fluxes of every level
    for (int lev = 0; lev < nlev; lev++){
        dUdt[lev].setVal(Real(0.0));
        FillPatchStageLevels(Uborder[lev], lev, stage_time, U_stage, 0, 0, NSTATE);

        if (stage == 0){
            if (Grad_Q){
                PhaseField_PreTimeStep(Uborder[lev], ncycle, stage_time, lev);
            }
            else {
                PreTimeStep(Uborder[lev], ncycle, stage_time, lev);
            }
        }

        DefineFaceFluxes(fluxes[lev], fluxes_nc[lev], lev);
        if (Grad_Q){
            PhaseField_ComputeFaceFluxes(Uborder[lev], (*Grad_Q)[lev], fluxes[lev], fluxes_nc[lev], lev, stage);
        }
        else {
            ComputeFaceFluxes(Uborder[lev], fluxes[lev], fluxes_nc[lev], lev, stage_time, stage);
            CheckCFL(lev, dt_stage);
        }
    }

    // 3. the coarse faces covered by the finer level take its fluxes and NC face quantities,
    //    averaged with the face areas (finest level first)
    if (do_reflux){
        for (int lev = finest_level-1; lev >= 0; lev--){
            AverageDownFacesTo(fluxes, lev);
#if (NONCONSERVATIVE == true)
            AverageDownFacesTo(fluxes_nc, lev);
#endif
        }
    }

    // 4. right-hand side of every level from its (possibly averaged) faces
    for (int lev = 0; lev < nlev; lev++){
        if (Grad_Q){
            PhaseField_FluxDivergence(dUdt[lev], Uborder[lev], fluxes[lev], fluxes_nc[lev], lev);
        }
        else {
            FluxDivergence(dUdt[lev], Uborder[lev], fluxes[lev], fluxes_nc[lev], lev, stage_time, dt_stage);
        }
    }
}


// One time step dt_lev of all levels, stage by stage (see the top of this file). mf_new and mf_old
// are swapped first, so that mf_old holds U^n; Grad_Q == nullptr: the hyperbolic step, otherwise
// the Phase-Field step. The stage combinations, stage times and stage weights are those of
// ForwardEuler, SecondOrderSSPRK, ThirdOrderSSPRK and FourthOrderRK. After each stage but the last
// the levels are averaged down; after the last, timeStepNoSubcycling (PhaseField_TimeStepNoSubcycling)
// does it
void
Compressible_PhaseField::AdvanceAllLevelsStageCoupled (Vector<MultiFab>& mf_new,
                                                       Vector<MultiFab>& mf_old,
                                                       Vector<Array<MultiFab,AMREX_SPACEDIM>> const* Grad_Q,
                                                       Real time,
                                                       Real dt_lev,
                                                       int ncycle)
{
    BL_PROFILE("AdvanceAllLevelsStageCoupled()");

    const int nlev = finest_level + 1;
    const int num_grow = NGROW;
    const bool hyperbolic = (Grad_Q == nullptr);

    for (int lev = 0; lev < nlev; lev++){
        t_old[lev] = t_new[lev];
        t_new[lev] += dt_lev;
    }

    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    for (int lev = 0; lev < nlev; lev++){
        std::swap(mf_old[lev], mf_new[lev]);
    }

    Vector<MultiFab> dUdt(nlev);
    Vector<MultiFab> Uborder(nlev);
    for (int lev = 0; lev < nlev; lev++){
        dUdt[lev].define(grids[lev], dmap[lev], NSTATE, 0);
        Uborder[lev].define(grids[lev], dmap[lev], NSTATE, num_grow);
    }

    auto rhs = [&] (Vector<MultiFab>& U_stage, Real stage_time, Real dt_stage, int stage)
    {
        StageRHSAllLevels(dUdt, Uborder, U_stage, Grad_Q, stage_time, dt_stage, ncycle, stage);
    };

    // energy-bound check of the hyperbolic step on every level; true: restart the step
    auto bound_violated = [&] (Vector<MultiFab>& U) -> bool
    {
        if (!hyperbolic) { return false; }
        for (int lev = 0; lev < nlev; lev++){
            if (StageBoundViolated(U[lev])) { return true; }
        }
        return false;
    };

    // end of a stage on every level (the relaxation of the six-equation models)
    auto post_stage = [&] (Vector<MultiFab>& U, Real t)
    {
        for (int lev = 0; lev < nlev; lev++){
            if (hyperbolic) {
                PostTimeStage(U[lev], ncycle, t, lev);
            }
            else {
                PhaseField_PostTimeStage(U[lev], ncycle, t, lev);
            }
        }
    };

    auto post_step = [&] (Vector<MultiFab>& U, Real t)
    {
        for (int lev = 0; lev < nlev; lev++){
            if (hyperbolic) {
                PostTimeStep(U[lev], ncycle, t, lev);
            }
            else {
                PhaseField_PostTimeStep(U[lev], ncycle, t, lev);
            }
        }
    };

    if (ID_TimeIntegrator == "ForwardEuler")
    {
        rhs(mf_old, time, dt_lev, 0);
        for (int lev = 0; lev < nlev; lev++){
            // U^* = U^n + dt*dUdt^n
            MultiFab::LinComb(mf_new[lev], Real(1.0), Uborder[lev], 0, dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time);
        post_step(mf_new, time);
    }
    else if (ID_TimeIntegrator == "TVD-RK2")
    {
        // RK2 stage 1
        rhs(mf_old, time, Real(0.5)*dt_lev, 0);
        for (int lev = 0; lev < nlev; lev++){
            // U^* = U^n + dt*dUdt^n
            MultiFab::LinComb(mf_new[lev], Real(1.0), Uborder[lev], 0, dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+dt_lev);
        AverageDown(mf_new);

        // RK2 stage 2, from U^* at time+dt
        rhs(mf_new, time+dt_lev, Real(0.5)*dt_lev, 1);
        for (int lev = 0; lev < nlev; lev++){
            // U_new = 0.5*(Uborder+U_old) + 0.5*dt*dUdt^*
            MultiFab::LinComb(mf_new[lev], Real(0.5), Uborder[lev], 0, Real(0.5), mf_old[lev], 0, 0, NSTATE, 0);
            MultiFab::Saxpy(mf_new[lev], Real(0.5)*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+dt_lev);
        post_step(mf_new, time+dt_lev);
    }
    else if (ID_TimeIntegrator == "TVD-RK3")
    {
        // RK3 stage 1
        rhs(mf_old, time, (Real(1.0)/Real(6.0))*dt_lev, 0);
        for (int lev = 0; lev < nlev; lev++){
            MultiFab::LinComb(mf_new[lev], Real(1.0), Uborder[lev], 0, dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+dt_lev);
        AverageDown(mf_new);

        // RK3 stage 2, from U^(1) at time+dt
        rhs(mf_new, time+dt_lev, (Real(1.0)/Real(6.0))*dt_lev, 1);
        for (int lev = 0; lev < nlev; lev++){
            MultiFab::LinComb(mf_new[lev], Real(0.25), Uborder[lev], 0, Real(0.75), mf_old[lev], 0, 0, NSTATE, 0);
            MultiFab::Saxpy(mf_new[lev], Real(0.25)*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+Real(0.5)*dt_lev);
        AverageDown(mf_new);

        // RK3 stage 3, from U^(2) at time+dt/2
        rhs(mf_new, time+Real(0.5)*dt_lev, (Real(2.0)/Real(3.0))*dt_lev, 2);
        for (int lev = 0; lev < nlev; lev++){
            MultiFab::LinComb(mf_new[lev], (Real(2.0)/Real(3.0)), Uborder[lev], 0, (Real(1.0)/Real(3.0)), mf_old[lev], 0, 0, NSTATE, 0);
            MultiFab::Saxpy(mf_new[lev], (Real(2.0)/Real(3.0))*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+dt_lev);
        post_step(mf_new, time+dt_lev);
    }
    else if (ID_TimeIntegrator == "RK4")
    {
        // the stage states u_1, u_2, u_3 are in U_temp, while mf_new holds the running sum
        Vector<MultiFab> U_temp(nlev);
        for (int lev = 0; lev < nlev; lev++){
            U_temp[lev].define(grids[lev], dmap[lev], NSTATE, 0);
        }

        // RK4 stage 1
        rhs(mf_old, time, (Real(1.0)/Real(6.0))*dt_lev, 0);
        for (int lev = 0; lev < nlev; lev++){
            /* u_1 = u_n + 1/2 dt R(u_n) */
            MultiFab::LinComb(U_temp[lev], Real(1.0), Uborder[lev], 0, Real(0.5)*dt_lev,             dUdt[lev], 0, 0, NSTATE, 0);
            /* u_* = u_n + 1/6 dt R(u_n) */
            MultiFab::LinComb(mf_new[lev], Real(1.0), Uborder[lev], 0, (Real(1.0)/Real(6.0))*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(U_temp)) { return; }
        post_stage(U_temp, time+Real(0.5)*dt_lev);
        AverageDown(U_temp);

        // RK4 stage 2, from u_1 at time+dt/2
        rhs(U_temp, time+Real(0.5)*dt_lev, (Real(1.0)/Real(3.0))*dt_lev, 1);
        for (int lev = 0; lev < nlev; lev++){
            /* u_2 = u_n + 1/2 dt R(u_1) */
            MultiFab::LinComb(U_temp[lev], Real(1.0), mf_old[lev], 0, Real(0.5)*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
            /* u_** = u_* + 1/3 dt R(u_1) */
            MultiFab::Saxpy(mf_new[lev], (Real(1.0)/Real(3.0))*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(U_temp)) { return; }
        post_stage(U_temp, time+Real(0.5)*dt_lev);
        AverageDown(U_temp);

        // RK4 stage 3, from u_2 at time+dt/2
        rhs(U_temp, time+Real(0.5)*dt_lev, (Real(1.0)/Real(3.0))*dt_lev, 2);
        for (int lev = 0; lev < nlev; lev++){
            /* u_3 = u_n + dt R(u_2) */
            MultiFab::LinComb(U_temp[lev], Real(1.0), mf_old[lev], 0, Real(1.0)*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
            /* u_*** = u_** + 1/3 dt R(u_2) */
            MultiFab::Saxpy(mf_new[lev], (Real(1.0)/Real(3.0))*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(U_temp)) { return; }
        post_stage(U_temp, time+dt_lev);
        AverageDown(U_temp);

        // RK4 stage 4, from u_3 at time+dt
        rhs(U_temp, time+dt_lev, (Real(1.0)/Real(6.0))*dt_lev, 3);
        for (int lev = 0; lev < nlev; lev++){
            /* u_n+1 = u_*** + 1/6 dt R(u_3) */
            MultiFab::Saxpy(mf_new[lev], (Real(1.0)/Real(6.0))*dt_lev, dUdt[lev], 0, 0, NSTATE, 0);
        }
        if (bound_violated(mf_new)) { return; }
        post_stage(mf_new, time+dt_lev);
        post_step(mf_new, time+dt_lev);
    }
}
