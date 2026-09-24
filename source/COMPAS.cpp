// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
// Portions derived from AMReX, BSD-3-Clause. See THIRD_PARTY_NOTICES.md.

#include <algorithm>
#include <numeric>
#include <vector>
#include <sstream>
#include <string>
#include <filesystem>
// #include <execution>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_TimeIntegrator.H>



namespace fs = std::filesystem;


#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <COMPAS.H>
#include <Kernels.H>
#include <Tools.H>
#include <IncludeOperators.H>

using namespace amrex;

// int FR_COMP = 2*NSTATE;
int FR_COMP = NSTATE;

Parm* COMPAS::h_parm = nullptr;
Parm* COMPAS::d_parm = nullptr;


// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
COMPAS::COMPAS ()
{

#ifdef PRECONSTRUCTOR
    PreConstructor();
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(vnamesC.size()) == NSTATE,
                                     "The state variable names must match the state components");
    
    variableSetUp();

    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    dof_new.resize(nlevs_max);
    dof_old.resize(nlevs_max);

    if (ID_IB){
        ib_levelset.resize(nlevs_max);
    }

    c_max.resize(nlevs_max);

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);
#if (NONCONSERVATIVE == true)
    flux_reg_nc.resize(nlevs_max+1);
#endif
}

COMPAS::~COMPAS ()
{
}


void
COMPAS::Physics_Derived_Constructor ()
{
    SetPhysicsBC();
}


void 
COMPAS::ExplicitTimeStep (amrex::Real cur_time)
{

    int lev = 0;
    int iteration = 1;
    int stage = 0;
    if (do_subcycle){
        timeStepWithSubcycling(dof_new, dof_old, lev, cur_time, iteration, stage);
    }
    else{
        timeStepNoSubcycling(cur_time, iteration);
    }
}

bool
COMPAS::OutputTimeInterval (amrex::Real cur_time)
{
    auto try_interval = [&](amrex::Real interval) -> bool
    {
        if (interval <= 0.0) return false;

        amrex::Real ratio_m =  cur_time          / interval;
        amrex::Real ratio_p = (cur_time + dt[0]) / interval;
        amrex::Real offset  = std::floor(ratio_p);

        ratio_m -= offset;
        ratio_p -= offset;

        amrex::Real dt_out = std::numeric_limits<amrex::Real>::infinity();

        if (ratio_m < -1.0e-12 && ratio_p > 0.0) {
            dt_out = offset * interval - cur_time;
        }

        if (interval < dt[0]) {
            dt_out = interval;
        }

        if (dt_out <= dt[0] && dt_out <= interval) {
            dt[0] = dt_out;
            for (int lev = 1; lev <= finest_level; ++lev) {
                dt[lev] = dt[lev-1] / nsubsteps[lev];
            }
            return true;
        }

        return false;
    };

    // 1. Try windowed output intervals first
    for (const auto& w : output_windows)
    {
        if (cur_time >= w.t_lo && cur_time < w.t_hi) {
            if (try_interval(w.dt)) return true;
        }
    }

    // 2. Fallback to global interval
    return try_interval(t_write_interval);
}


void
COMPAS::TimeStepChecks(int step)
{
    if (step >= cfl_switch){
        cfl = cfl_fast;
    }

    TimeStepChecksDerived(step);
}

void
COMPAS::SumVar()
{
            // sum a variable to check conservation and for NaNs 
    sum_state = dof_new[0].sum(sum_var);
    
    if (std::isnan(sum_state)){
#ifdef USER_OUTPUT_FUNC 
        UserOutputFunction(*this);
#endif
        WriteCheckpointFile();
        amrex::Abort("nan found in the state!\n");
    }
}

void 
COMPAS::TimestepOutput(int step, amrex::Real cur_time)
{
    amrex::Print() << "[COMPAS] Coarse step:     " << step << "\n";
    amrex::Print() << "           SimulationTime:  " << cur_time << "\n";
    amrex::Print() << "           Timestep:        " << dt[0] << "\n";

    if (timestep_sum){
    amrex::Print() << "           SumVar:          " << vnamesC[sum_var] << "\n";
    amrex::Print() << "           Sum:             " << sum_state << "\n";
    }

    if (timers){

        amrex::Real est_tot_time = evolve_total + (stop_time - cur_time)/dt[0]*evolve_step;

        amrex::Real ra_est_tot_time = 0.0;
        //ra_estimates.push_back(est_tot_time);
        //if (step <= run_time_ra_int){
        //    //ra_est_tot_time = std::reduce(ra_estimates.begin(), ra_estimates.end())/step;
        //    ra_est_tot_time =
	//	        std::reduce(ra_estimates.begin(),
	//				                ra_estimates.end(),
	//						                0.0) / step;
	//} else {
        //    ra_estimates.erase(ra_estimates.begin());
            //ra_est_tot_time = std::reduce(ra_estimates.begin(), ra_estimates.end())/run_time_ra_int;
       	//    ra_est_tot_time =
	//	        std::reduce(ra_estimates.begin(),
	//				                ra_estimates.end(),
	//						                0.0) / ra_time_ra_int;
	//}
        // amrex::Print() << "step: " << step << " ra_est_tot_time:" << ra_est_tot_time << "\n";
        // for (int i = 0; i < ra_estimates.size(); i++){
        //     amrex::Print() << "ra_estimates[" << i << "]: " << ra_estimates[i] << "\n";
        // }
        
        std::string tot_time_units = "seconds";
        if (est_tot_time > 60.0 && est_tot_time < 3600.0){
            est_tot_time /= 60.0;
            tot_time_units = "minutes";
        } else if (est_tot_time > 3600.0){
            est_tot_time /= 3600.0;
            tot_time_units = "hours";
        }

        std::string ra_tot_time_units = "seconds";
        if (ra_est_tot_time > 60.0 && ra_est_tot_time < 3600.0){
            ra_est_tot_time /= 60.0;
            ra_tot_time_units = "minutes";
        } else if (ra_est_tot_time > 3600.0){
            ra_est_tot_time /= 3600.0;
            ra_tot_time_units = "hours";
        }

        amrex::Real cur_run_time = evolve_total;
        std::string cur_time_units = "seconds";
        if (cur_run_time > 60.0 && cur_run_time < 3600.0){
            cur_run_time /= 60.0;
            cur_time_units = "minutes";
        } else if (cur_run_time > 3600.0){
            cur_run_time /= 3600.0;
            cur_time_units = "hours";
        }

        std::string ra_int_order_str = std::to_string(run_time_ra_int);
        int ra_int_order = std::size(ra_int_order_str);


        amrex::Print() << "           EvolveTime:    " << std::string(ra_int_order+3, ' ') << evolve_step << " seconds\n";
        amrex::Print() << "           EstTotalTime:  " << std::string(ra_int_order+3, ' ') << est_tot_time << " " << tot_time_units << "\n";
        amrex::Print() << "           EstTotalTimeRA(" << run_time_ra_int << "): " << ra_est_tot_time << " " << ra_tot_time_units << "\n";
        amrex::Print() << "           CurrentRunTime:" << std::string(ra_int_order+3, ' ') << cur_run_time << " " << cur_time_units << "\n";
    }
            
}


// advance solution to final time
void
COMPAS::Evolve ()
{

#if (CONVERGENCE == true)
    Evolve_start_time = amrex::second();
#endif

    Real cur_time = t_new[0];
    int last_plot_file_step = 0;
    int twi = 0;
    int lev = 0;

    amrex::Real evolve_start = amrex::second();
    if (timers){
        evolve_total = 0.0;
    }

    /* Loop over time steps, stopping when the final time is reached */
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\n==== Coarse step " << step+1 << " ====" << std::endl;

// #if (OUTPUT_TOOLS == true)
//         CalculateOutputMF();
//         if ((step+1) % save_points_int == 0 && save_points_int > 0){
//             get_data_at_points (dof_new);
//         }
//         if ((step+1) % save_lines_int == 0 && save_lines_int > 0){
//             get_data_at_lines (dof_new);
//         }
//         ResetOutputMF();
// #endif

        TimeStepChecks(step);

        /* Compute the time step on each level */
#if (FIXED_DT == true)
        FixedDt();
#else
        ComputeDt();
#endif

        bool OutputInterval = OutputTimeInterval(cur_time);

        if (timers){
            evolve_step = amrex::second();
        }

        ExplicitTimeStep(cur_time);

        cur_time += dt[0];

        if (timers){
            evolve_step = amrex::second() - evolve_step;
            evolve_total = amrex::second() - evolve_start;
        }


        //=================================================================================
        //                           Post timestep things
        //=================================================================================
        if (timestep_sum){
            SumVar();
        }

        TimestepOutput(step+1,cur_time);

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

#ifdef USER_OUTPUT_FUNC 
        if (user_output_int > 0 && (step+1) % user_output_int == 0) {
            UserOutputFunction(*this);
        }
#endif

        // One plotfile for the listed times this step reached. Every listed time already reached
        // is passed over, so that two times in one step or times before a restart do not stop the rest
        bool OutputTime = false;
        while (twi < t_write.size() && t_write[twi] <= cur_time) {
            OutputTime = OutputTime || t_write[twi] >= cur_time - dt[0];
            twi++;
        }

        if (
            (plot_int > 0 && (step+1) % plot_int == 0) || 
            OutputTime
            || OutputInterval
            ) {
            last_plot_file_step = step+1;

            CalculateOutputMF();
            WritePlotFile();
            ResetOutputMF();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (timers){
        avg_timestep = evolve_total/(istep[0] + 1);    
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {

#ifdef USER_OUTPUT_FUNC 
        if (user_output_int <= 0) {
            UserOutputFunction(*this);
        }
#endif

        CalculateOutputMF();
        WritePlotFile();
// #if (OUTPUT_TOOLS == true)
//         if (save_points_int > 0){
//             get_data_at_points (dof_new);
//         }
//         if (save_lines_int > 0){
//             get_data_at_lines (dof_new);
//         }
// #endif
        ResetOutputMF();

        if (chk_int > 0 ) {
            WriteCheckpointFile();
        }
    }
    
    /* This was causing some weird behavior - investigate */
    // variableCleanUp();

#if (CONVERGENCE == true)
    Evolve_stop_time = amrex::second() - Evolve_start_time;
    OutputError();
#endif

}

// initializes multilevel data
void
COMPAS::InitData ()
{

    BL_PROFILE("InitData()");

    Physics_Derived_Constructor();

    /* IT IS VERY IMPORTANT THAT DOF_NEW HAS THE PROPER DATA HERE */
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown(dof_new);

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
        
        if (istep[0] >= cfl_switch){
            cfl = cfl_fast;
        }
        // Get c_max for the time-step calculation
        for (int lev = 0; lev <= finest_level; lev++){
            const Real* dx  =  geom[lev].CellSize();

#if (ADVECTION == true)            
            amrex::Real c_max_est_adv = 0.0;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                c_max_est_adv = std::max(c_max_est_adv,cfl*dx[idim]/dt[0]);
            }
            c_max[lev].setVal(c_max_est_adv,0,1,0);
#endif
#if (DIFFUSION==true)    
            amrex::Real c_max_est_diff = 0.0;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                c_max_est_diff = std::max(c_max_est_diff,vnn*dx[idim]*dx[idim]/dt[0]);
            }
            c_max[lev].setVal(c_max_est_diff,1,1,0);        
#endif
        }
    }

    if (plot_int > 0) {
        CalculateOutputMF();
        WritePlotFile();
        ResetOutputMF();
    }
}


// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
COMPAS::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{

    BL_PROFILE("MakeNewLevelFromCoarse()");

    /* IT IS VERY IMPORTANT THAT DOF_NEW HAS THE PROPER DATA HERE */
    const int ncomp = NSTATE;
    const int ncomp_nc = NC_TERMS; // Change if needed
    const int nghost = dof_new[lev-1].nGrow();

    dof_new[lev].define(ba, dm, ncomp, nghost);
    dof_old[lev].define(ba, dm, ncomp, nghost);

    if (ID_IB){
        ib_levelset[lev].define(ba, dm, 1, nghost);
    }
    

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    //define the maximum speed with the box array, distribution mapping, one component, and zero
    int n_cmax_comp = (ADVECTION == true) + (DIFFUSION == true);
    c_max[lev].define(ba, dm, n_cmax_comp, 1);

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, FR_COMP));
#if (NONCONSERVATIVE == true)
        flux_reg_nc[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_nc));
#endif
    }

    /* CHANGE MF HERE */
    FillCoarsePatch(dof_new[lev], lev, time,
                    dof_new, t_new,
                    dof_old, t_old,
                    0, 0, ncomp);

    MakeNewLevelFromCoarse_Derived(lev,time,ba,dm);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
COMPAS::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{

    BL_PROFILE("RemakeLevel()");

    /* IT IS VERY IMPORTANT THAT DOF_NEW HAS THE PROPER DATA HERE */
    const int ncomp = NSTATE;
    const int ncomp_nc = NC_TERMS; // Change if needed
    const int nghost = dof_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    

    FillPatch(new_state, lev, time,
              dof_new, t_new,
              dof_old, t_old,
              0, 0, ncomp);

    std::swap(new_state, dof_new[lev]);
    std::swap(old_state, dof_old[lev]);


    if (ID_IB){
        MultiFab new_ib_levelset(ba, dm, 1, nghost);
        FillPatch(new_ib_levelset, lev, time,
              ib_levelset, t_new,
              ib_levelset, t_old,
              0, 0, 1);
        std::swap(new_ib_levelset, ib_levelset[lev]);
    }


    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    // for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    // {
    //     facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
    // }

    int n_cmax_comp = (ADVECTION == true) + (DIFFUSION == true);
    c_max[lev] = MultiFab(ba, dm, n_cmax_comp, 1);

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, FR_COMP));
#if (NONCONSERVATIVE == true)
        flux_reg_nc[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_nc));
#endif
    }

    RemakeLevel_Derived(lev,time,ba,dm);
}


// Remake all levels using provided input data
void
COMPAS::RemakeAllLevelsFromData (Vector<MultiFab> const& dof_Input,
                                   Vector<Real>     const&   t_Input,
                                   Vector<Real>     const&  dt_Input,
                                   Vector<int>      const& istep_Input,
                                   int              finest_level_Input)
{
    //Reset finest_level
    finest_level = finest_level_Input;

    //Reset istep to be all zero
    istep.assign(max_level+1, 0);

    for (int lev = 0; lev <= finest_level; lev++){
        //Reset BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray       (lev, dof_Input[lev].boxArray());
        SetDistributionMap(lev, dof_Input[lev].DistributionMap());

        //Reset dof_new and dof_old
        int Num_Comp = dof_Input[lev].nComp();
        int Num_Grow = dof_Input[lev].nGrow();
        dof_new[lev].define(grids[lev], dmap[lev], Num_Comp, Num_Grow);
        dof_old[lev].define(grids[lev], dmap[lev], Num_Comp, Num_Grow);
        MultiFab::Copy(dof_new[lev], dof_Input[lev], 0, 0, Num_Comp, Num_Grow);
        MultiFab::Copy(dof_old[lev], dof_Input[lev], 0, 0, Num_Comp, Num_Grow);

        //Reset FluxRegister
        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, FR_COMP));
#if (NONCONSERVATIVE == true)
            flux_reg_nc[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, NC_TERMS));
#endif
        }

        //Reset c_max
        Num_Comp = (ADVECTION == true) + (DIFFUSION == true);
        c_max[lev].define(grids[lev], dmap[lev], Num_Comp, 1);
        c_max[lev].setVal(0.0);

        //Reset t_new
        t_new[lev] = t_Input[lev];

        //Reset dt
        dt[lev] = dt_Input[lev];

        //Reset istep
        istep[lev] = istep_Input[lev];

        //Reset model-specific data
        RemakeAllLevelsFromData_Derived(lev, grids[lev], dmap[lev]);
    }
    
}


// Delete level data
// overrides the pure virtual function in AmrCore
void
COMPAS::ClearLevel (int lev)
{
    dof_new[lev].clear();
    dof_old[lev].clear();
    c_max[lev].clear();
    flux_reg[lev].reset(nullptr);
    if (ID_IB){
        ib_levelset[lev].clear();
    }
#if (NONCONSERVATIVE == true)
    flux_reg_nc[lev].reset(nullptr);
#endif

    ClearLevel_Derived(lev);
}

amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> COMPAS::dxFinest()
{
    auto dx_Finest = Geom(0).CellSizeArray(); 
    for (int lev = 1; lev < max_level + 1; lev++){
        for (int dim = 0; dim < AMREX_SPACEDIM; dim++){
            auto rr = refRatio(lev-1);
            dx_Finest[dim] /= rr[dim];
        }
    }
    return dx_Finest;
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void COMPAS::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{

    BL_PROFILE("MakeNewLevelFromScratch()");

    /* Assume this function does have anything to do with time
      stepping, and dof_new can be filled as such */

    /* IT IS VERY IMPORTANT THAT DOF_NEW HAS THE PROPER DATA HERE */

    const int ncomp = NSTATE;
    const int ncomp_nc = NC_TERMS; // Change if needed

    dof_new[lev].define(ba, dm, ncomp, 1);
    dof_old[lev].define(ba, dm, ncomp, 1);

    if (ID_IB){
        ib_levelset[lev].define(ba, dm, 1, 1);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
        flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, FR_COMP));
#if (NONCONSERVATIVE == true)
        flux_reg_nc[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp_nc));
#endif
    }

    int n_cmax_comp = (ADVECTION == true) + (DIFFUSION == true);
    c_max[lev].define(ba, dm, n_cmax_comp, 1);

    MakeNewLevelFromScratch_Derived(lev,time,ba,dm);

    MultiFab& state_old = dof_old[lev];
    MultiFab& state_new = dof_new[lev];

    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    const auto dx_FinestLevel = dxFinest();

    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state_old,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> fabo = state_old[mfi].array();
        Array4<Real> fabn = state_new[mfi].array();
        const Box& box = mfi.tilebox();

        amrex::launch(box,
        [=] AMREX_GPU_DEVICE (Box const& tbx)
        {
            initdata(tbx, fabn, problo, dx, dx_FinestLevel, ncomp, *lparm, ID_QuadratureIC);
        });

        if (ID_IB){
            Array4<Real> fab_ib = ib_levelset[lev][mfi].array();
            amrex::launch(box,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                initdata(tbx, fab_ib, problo, dx, dx_FinestLevel, 1, *lparm, ID_QuadratureIC, 1);
            });
        }

        fabo = fabn;
    }

}


// overrides the pure virtual function in AmrCore
void
COMPAS::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    BL_PROFILE("ErrorEst()");

    MultiFab& state = dof_new[lev];

    // State with ghost cells. Derived variables are computed on ngrow ghost cells, and the
    // ones from the velocity gradient need one more
    int ngrow = 1;
    int ngrow_state = ngrow;
    for (const TagVar& TV : tag_vars) {
        if (TV.derived && DerivedVarNeedsNeighbors(TV.name)) ngrow_state = ngrow + 1;
    }
    MultiFab Sborder(grids[lev], dmap[lev], NSTATE, ngrow_state);
    FillPatch(Sborder, lev, t_new[lev],
              dof_new, t_new,
              dof_old, t_old,
              0, 0, NSTATE);

    // Scratch space for derived tagging variables, kept apart from the state so that
    // computing one derived variable does not overwrite the state used by the next
    MultiFab Derived(grids[lev], dmap[lev], 1, ngrow);

    const Real dx = geom[lev].CellSize(0);
#if (AMREX_SPACEDIM > 1)
    const Real dy = geom[lev].CellSize(1);
#else
    const Real dy = 1.0;
#endif
#if (AMREX_SPACEDIM > 2)
    const Real dz = geom[lev].CellSize(2);
#else
    const Real dz = 1.0;
#endif

    const Real time = t_new[lev];
    const auto problo = Geom(lev).ProbLoArray();
    const auto dX     = Geom(lev).CellSizeArray();

    const int tagval = TagBox::SET;
    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const Box& gbx = amrex::grow(bx, ngrow);

            const auto statefab   = Sborder.array(mfi);
            const auto derivedfab = Derived.array(mfi);
            const auto tagfab     = tags.array(mfi);

            // Array tagged below, either the state or the derived variable
            Array4<Real> varfab = statefab;

            if (ID_IB){
                const auto ib_fab   = ib_levelset[lev].array(mfi);
                Real val_low = -0.1;
                Real val_high = 0.1;
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    state_value_range(i, j, k, 0, tagfab, ib_fab, val_low, val_high, tagval);
                });
            }
            
            // Loop over refinement variables
            for (int v = 0; v < tag_vars.size(); ++v)
            {
                const TagVar& TV = tag_vars[v];

                // Enforce max refinement level
                if (lev >= TV.max_level) continue;

                int comp = TV.index;

                // --- Derived variable handling ---
                if (TV.derived)
                {
                    // Compute derived quantity into component 0
                    amrex::ParallelFor(gbx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        compute_derived_var(i, j, k,
                                             comp,
                                             dx, dy, dz,
                                             derivedfab,
                                             statefab,
                                             *lparm);
                    });

                    varfab = derivedfab;
                    comp = 0;  // derived variable always tagged from comp 0
                }
                else
                {
                    // Use state data directly
                    varfab = statefab;
                }

                // --- Gradient-based tagging ---
                if (TV.do_grad)
                {
                    Real val = TV.grad[lev];
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        state_grad(i, j, k, comp, tagfab, varfab, val, tagval);
                    });
                }

                // --- Value-less-than tagging ---
                if (TV.do_value_less)
                {
                    Real val = TV.value_less[lev];
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        state_value_less(i, j, k, comp, tagfab, varfab, val, tagval);
                    });
                }

                // --- Value-greater-than tagging ---
                if (TV.do_value_greater)
                {
                    Real val = TV.value_greater[lev];
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        state_value_greater(i, j, k, comp, tagfab, varfab, val, tagval);
                    });
                }

                // --- Value-range tagging ---
                if (TV.do_value_range)
                {
                    Real val_low  = TV.value_range[2*lev];
                    Real val_high = TV.value_range[2*lev + 1];
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        state_value_range(i, j, k,
                                          comp, tagfab, varfab,
                                          val_low, val_high,
                                          tagval);
                    });
                }
            }

#ifdef SPACE_TIME_REF
            // Spatio-temporal tagging, on the valid cells like the criteria above, because
            // the tags have no ghost cells when amr.n_error_buf = 0
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                space_time(i, j, k, lev, time, problo, dX, tagfab, *lparm);
            });
#endif
        }
    }
}

void 
COMPAS::ResolvePaths()
{
    auto split_underscore = [](const std::string& s) {
        std::vector<std::string> parts;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, '_')) {
            if (!item.empty()) parts.push_back(item);
        }
        return parts;
    };

    auto dirname = [](const std::string& path) {
        size_t pos = path.find_last_of('/');
        if (pos == std::string::npos) return std::string{};
        if (pos == 0) return std::string{"/"};
        return path.substr(0, pos);
    };

    auto basename = [](const std::string& path) {
        size_t pos = path.find_last_of('/');
        if (pos == std::string::npos) return path;
        return path.substr(pos + 1);
    };

    auto join_path = [](const std::string& a, const std::string& b) {
        if (a.empty() || a == ".") return b;
        if (a.back() == '/') return a + b;
        return a + "/" + b;
    };

    // ======================================================
    // 1. Safety check
    // ======================================================
    if (!output_dir.empty() && !case_name.empty() && !plot_file.empty()) {
        amrex::Abort(
            "Error: Cannot specify output_dir, case_name, and plot_file simultaneously.\n"
        );
    }

    // ======================================================
    // 2. User provided ONLY case_name (may include dirs)
    // ======================================================
    if (!case_name.empty() && output_dir.empty() && plot_file.empty()) {

        std::string cname = basename(case_name);
        if (cname.empty()) {
            amrex::Abort("Error: case_name path is malformed.\n");
        }

        std::string parent = dirname(case_name);
        output_dir = parent.empty() ? "." : parent;

        case_name = cname;

        std::string case_dir = join_path(output_dir, case_name);
        plot_file      = join_path(case_dir, "plt_" + case_name);
        movie_filename = join_path(case_dir, case_name + ".visit");

        return;
    }

    // ======================================================
    // 3. User provided output_dir + case_name
    // ======================================================
    if (!output_dir.empty() && !case_name.empty() && plot_file.empty()) {

        std::string case_dir = join_path(output_dir, case_name);
        plot_file      = join_path(case_dir, "plt_" + case_name);
        movie_filename = join_path(case_dir, case_name + ".visit");

        return;
    }

    // ======================================================
    // 4. Otherwise parse plot_file
    // ======================================================
    if (plot_file.empty()) {
        amrex::Abort(
            "Error: plot_file must be provided when case_name/output_dir are not.\n"
        );
    }

    output_dir = dirname(plot_file);
    if (output_dir.empty()) output_dir = ".";

    std::string filename = basename(plot_file);

    auto tokens = split_underscore(filename);

    if (tokens.size() < 2) {
        amrex::Abort(
            "Error: plot_file must look like: plt_CASE or plt_CASE_XXXXX.\n"
        );
    }

    if (tokens[0].substr(0, 3) != "plt") {
        amrex::Abort("Error: plot_file must begin with 'plt'.\n");
    }

    bool last_is_number = std::all_of(
        tokens.back().begin(), tokens.back().end(), ::isdigit
    );

    size_t end_idx = last_is_number ? tokens.size() - 1 : tokens.size();

    case_name.clear();
    for (size_t i = 1; i < end_idx; ++i) {
        case_name += tokens[i];
        if (i + 1 < end_idx) case_name += "_";
    }

    if (case_name.empty()) {
        amrex::Abort(
            "Error: Could not extract case_name from plot_file.\n"
        );
    }

    // ======================================================
    // 5. Final output names
    // ======================================================
    {
        std::string case_dir = join_path(output_dir, case_name);
        plot_file      = join_path(case_dir, "plt_" + case_name);
        movie_filename = join_path(case_dir, case_name + ".visit");
    }
}


// void 
// COMPAS::ResolvePaths()
// {
//     auto split_underscore = [](const std::string& s) {
//         std::vector<std::string> parts;
//         std::stringstream ss(s);
//         std::string item;
//         while (std::getline(ss, item, '_')) {
//             if (!item.empty()) parts.push_back(item);
//         }
//         return parts;
//     };

//     // ======================================================
//     // 1. Safety check: user provided all 3.
//     // ======================================================
//     if (!output_dir.empty() && !case_name.empty() && !plot_file.empty()) {
//         amrex::Abort("Error: Cannot specify output_dir, case_name, and plot_file simultaneously.\n");
//     }

//     // ======================================================
//     // 2. User provided ONLY case_name
//     //    (New: case_name may include directories!)
//     // ======================================================
//     if (!case_name.empty() && output_dir.empty() && plot_file.empty()) {

//         fs::path case_path(case_name);

//         // Extract base case name (e.g. "advection_test")
//         std::string cname = case_path.filename().string();
//         if (cname.empty()) {
//             amrex::Abort("Error: case_name path is malformed.\n");
//         }

//         // Parent directory becomes output_dir (or "." if none)
//         fs::path parent = case_path.parent_path();
//         output_dir = parent.empty() ? "." : parent.string();

//         // Normalize stored case_name
//         case_name = cname;

//         // Construct final directory & filenames
//         fs::path case_dir = fs::path(output_dir) / case_name;
//         plot_file      = (case_dir / ("plt_" + case_name)).string();
//         movie_filename = (case_dir / (case_name + ".visit")).string();

//         return;
//     }

//     // ======================================================
//     // 3. User provided output_dir + case_name (no slashes).
//     // ======================================================
//     if (!output_dir.empty() && !case_name.empty() && plot_file.empty()) {

//         fs::path case_dir = fs::path(output_dir) / case_name;

//         plot_file      = (case_dir / ("plt_" + case_name)).string();
//         movie_filename = (case_dir / (case_name + ".visit")).string();

//         return;
//     }

//     // ======================================================
//     // 4. Otherwise parse plot_file
//     // ======================================================
//     if (plot_file.empty()) {
//         amrex::Abort(
//             "Error: plot_file must be provided when case_name/output_dir are not.\n"
//         );
//     }

//     fs::path p(plot_file);

//     // Extract directory
//     fs::path dir = p.parent_path();
//     output_dir = dir.empty() ? "." : dir.string();

//     // Filename (e.g. plt_case_00010)
//     std::string filename = p.filename().string();

//     // Split by underscores
//     auto tokens = split_underscore(filename);

//     if (tokens.size() < 2) {
//         amrex::Abort("Error: plot_file must look like: plt_CASE or plt_CASE_XXXXX.\n");
//     }

//     if (tokens[0].substr(0,3) != "plt") {
//         amrex::Abort("Error: plot_file must begin with 'plt'.\n");
//     }

//     // Check if last token is a timestep index
//     bool last_is_number = std::all_of(
//         tokens.back().begin(), tokens.back().end(), ::isdigit
//     );

//     size_t end_idx = last_is_number ? tokens.size() - 1 : tokens.size();

//     // Build case_name = tokens[1..end_idx-1], joined with underscores
//     case_name.clear();
//     for (size_t i = 1; i < end_idx; ++i) {
//         case_name += tokens[i];
//         if (i + 1 < end_idx) case_name += "_";
//     }

//     if (case_name.empty()) {
//         amrex::Abort("Error: Could not extract case_name from plot_file.\n");
//     }

//     // ======================================================
//     // 5. Final output names (unified formatting)
//     // ======================================================
//     {
//         fs::path case_dir = fs::path(output_dir) / case_name;

//         plot_file      = (case_dir / ("plt_" + case_name)).string();
//         movie_filename = (case_dir / (case_name + ".visit")).string();
//     }
// }

// read in some parameters from inputs file
void
COMPAS::ReadParameters ()
{
    {
        ParmParse pp("convergence");
        pp.query("dt_fixed",dt_fixed);
        pp.query("conv_output_file",conv_output_file);
        pp.query("conv_prob_name",conv_prob_name);

    }

    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("output_dir", output_dir);
        pp.query("case_name", case_name);
        pp.query("output_subdirectory", output_subdirectory);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
        pp.queryarr("t_write",t_write);
        pp.query("t_write_interval",t_write_interval);
        pp.queryarr("n_cell",n_cells);
        pp.query("Interpolater",ID_Interpolater);
        if (ID_Interpolater != "pc_interp" && ID_Interpolater != "lincc_interp" && ID_Interpolater != "quartic_interp"){
            amrex::Abort("amr.Interpolater = '" + ID_Interpolater + "' is not recognized. "
                         "Valid values are pc_interp, lincc_interp, quartic_interp.");
        }


        // Global fallback interval
        pp.query("t_write_interval", t_write_interval);

        int nwin = 0;
        pp.query("n_output_windows", nwin);

        output_windows.resize(nwin);

        for (int i = 0; i < nwin; ++i)
        {
            std::string prefix = "amr.output_window_" + std::to_string(i);
            amrex::ParmParse ppw(prefix);

            ppw.get("t_lo", output_windows[i].t_lo);
            ppw.get("t_hi", output_windows[i].t_hi);
            ppw.get("dt",   output_windows[i].dt);
        }

        if (t_write.size() == 0){
            t_write.resize(1); t_write[0] = -1.0;
        }
    }

    {
        ResolvePaths();
        
        std::vector<std::pair<std::string,std::string>> files = {
            {"ProblemICBC.H", "./prob/ProblemICBC.H"},
            {"inputs",        "./prob/inputs"},
            {"Parm.H",        "./prob/Parm.H"},
            {"UserParm.cpp",  "./prob/UserParm.cpp"},
            {"GNUmakefile",   "./GNUmakefile"}
        };

        DumpConfigurationLog(
            this->plot_file,   // from ResolvePaths()
            this->case_name,   // from ResolvePaths()
            files              // generalized list of sections
        );
    }

    {
        ParmParse pp("run");

        pp.query("cfl", cfl);
        pp.query("cfl_fast", cfl_fast);
        pp.query("cfl_switch", cfl_switch);
        pp.query("timestep_change_limiter",timestep_change_limiter);
        pp.query("vnn", vnn);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);
        pp.queryarr("output_vars", output_vars);
        pp.queryarr("lo_bc", lo_bc);
        pp.queryarr("hi_bc", hi_bc);
        if (lo_bc.size() < AMREX_SPACEDIM || hi_bc.size() < AMREX_SPACEDIM){
            amrex::Abort("run.lo_bc and run.hi_bc need one value per direction.");
        }
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            for (int bc : {lo_bc[idim], hi_bc[idim]}){
                if (bc < 0 || bc > 5){
                    amrex::Abort("run.lo_bc/run.hi_bc = " + std::to_string(bc) + " is not recognized. Valid values are "
                                 "0 (interior/periodic), 1 (inflow), 2 (outflow), 3 (symmetry), 4 (slip wall), 5 (no-slip wall).");
                }
            }
            bool periodic = Geom(0).isPeriodic(idim);
            if (periodic != (lo_bc[idim] == 0) || periodic != (hi_bc[idim] == 0)){
                amrex::Abort("Direction " + std::to_string(idim) + ": run.lo_bc/run.hi_bc = 0 must be used exactly when "
                             "geometry.is_periodic = 1.");
            }
        }
        pp.query("check_cfl",check_cfl);

        pp.query("QuadratureIC", ID_QuadratureIC);

        pp.query("TimeIntegrator", ID_TimeIntegrator);
        if (ID_TimeIntegrator != "ForwardEuler" && ID_TimeIntegrator != "TVD-RK2" &&
            ID_TimeIntegrator != "TVD-RK3" && ID_TimeIntegrator != "RK4"){
            amrex::Abort("run.TimeIntegrator = '" + ID_TimeIntegrator + "' is not recognized. "
                         "Valid values are ForwardEuler, TVD-RK2, TVD-RK3, RK4.");
        }
        pp.query("timestep_sum",timestep_sum);
        pp.query("user_output_int", user_output_int);
        pp.query("timers",timers);
        pp.query("sum_var", sum_var);
        pp.query("run_time_ra_int",run_time_ra_int);

        

        if (!cfl_fast){
            cfl_fast = cfl;
        }
        if (!cfl_switch){
            cfl_switch = 10;
        }

        if (output_vars.size() > 0){
            output_state_vars = false;
            num_output_vars = output_vars.size();
            if (num_output_vars > 2*NSTATE){
                amrex::Abort("Number of output vars is too large! If you would like to output more than 2*NSTATE variables, array size needs to change!\n");
            }

            for (int v = 0; v < num_output_vars; v++){
                std::vector<std::string>::iterator itr = std::find(vnamesD.begin(), vnamesD.end(), output_vars[v]);
                if (itr != vnamesD.cend()) {
                    out_vars[v] = std::distance(vnamesD.begin(), itr);
                    out_vars_string.push_back(vnamesD[out_vars[v]]);
                } else {
                    out_vars[v] = -1;
                    amrex::Print() << "I failed to find the variable " << output_vars[v] << "\n";
                    amrex::Abort("");
                }
            }
        } else {
            for (int v = 0; v < NSTATE; v++){
                out_vars[v] = v;
            }
            out_vars_string = vnamesC;
        }

    }

    {
        ParmParse pp("geometry");

        pp.query("coord_sys", ID_COORDSYS);
        if (ID_COORDSYS != 0){
            amrex::Abort("geometry.coord_sys = " + std::to_string(ID_COORDSYS) + " is not supported. "
                         "COMPAS supports Cartesian geometry only (geometry.coord_sys = 0).");
        }
    }

    
    {
        int proc = amrex::ParallelDescriptor::MyProc();
        
        ParmParse pp("save.lines");
        pp.queryarr("z_lines_to_save",z_lines_to_save);
        pp.query("save_lines_int",save_lines_int);
        pp.query("filename",lines_filename);

        num_save_straight_lines = int(z_lines_to_save.size()/2);

        /* Clean out the files for each save point */
        for (int p = 0; p < num_save_straight_lines; p++){
            Real px = z_lines_to_save[0 + 2*p];
            Real py = z_lines_to_save[1 + 2*p];

            for (int v = 0; v < num_output_vars; v++){
                {std::ofstream fp;
                fp.open (lines_filename + 
                        AMREX_D_TERM("_" + std::to_string(px) +,
                                    "_" + std::to_string(py) +,
                                    "_" + vnamesD[out_vars[v]] + "_" + std::to_string(proc) +) ".csv", std::ios::trunc);
                fp.precision(16);
                // fp << "t,";
                // for (int v = 0; v < num_output_vars; v++){
                // fp << vnamesD[out_vars[v]] << ",";

                // std::ofstream fp;
                //                     fp.open (lines_filename + 
                //                             AMREX_D_TERM("_" + std::to_string(px) +,
                //                                         "_" + std::to_string(py) +,
                //                                         "_" + vnamesD[out_vars[v]] + "_z" +) ".csv", std::ios::trunc);
                //                     fp.precision(16);
                                    
                }
            }
            // fp << "\n";
            {std::ofstream fp;
                fp.open (lines_filename + 
                        AMREX_D_TERM("_" + std::to_string(px) +,
                                    "_" + std::to_string(py) +,
                                    "_t" +) ".csv", std::ios::trunc);
                fp.precision(16);
                // fp << "t,";
            }
                {std::ofstream fp;
                                fp.open (lines_filename + 
                                        AMREX_D_TERM("_" + std::to_string(px) +,
                                                    "_" + std::to_string(py) +,
                                                     "_z_" + std::to_string(proc) +) ".csv", std::ios::trunc);
                                fp.precision(16);
                                // fp << "x" << ",";
                                // fp << "y" << ",";
                                // fp << "z" << ",";
                                // fp << "\n";
                                }
        }

        if (save_lines_int == 0){
            save_lines_int = 1000000000;
        }

    }


    {
        ParmParse pp("save.points");

        pp.queryarr("points_to_save",points_to_save);
        pp.query("save_points_int",save_points_int);
        pp.query("filename",points_filename);

        num_save_points = int(points_to_save.size()/3);

        /* Clean out the files for each save point */
        for (int p = 0; p < num_save_points; p++){
            Real px = points_to_save[0 + 3*p];
            Real py = points_to_save[1 + 3*p];
            Real pz = points_to_save[2 + 3*p];
            std::ofstream fp;
            fp.open (points_filename + 
                     AMREX_D_TERM("_" + std::to_string(px) +,
                                  "_" + std::to_string(py) +,
                                  "_" + std::to_string(pz) +) ".csv", std::ios::trunc);
            fp.precision(16);
            fp << "t,";
            for (int v = 0; v < num_output_vars; v++){
                fp << vnamesD[out_vars[v]] << ",";
            }
            fp << "\n";
        }

        if (save_points_int == 0){
            save_points_int = 1000000000;
        }

    }

    {
        int n = -1;
        ParmParse pp("run.IO");
        pp.query("n_out_files",n);
        if (n != -1){
            amrex::VisMF::SetNOutFiles(n);    
        }
        
    }

    {
        ParmParse pp("run.refine");

        amrex::Vector<std::string> ref_vars;
        pp.queryarr("ref_vars", ref_vars);

        tag_vars.resize(ref_vars.size());

        const auto& registry = VariableRegistry::instance();

        for (int i = 0; i < ref_vars.size(); ++i)
        {
            const std::string& name = ref_vars[i];

            if (!registry.has(name)) {
                amrex::Abort("Unknown refinement variable: " + name);
            }

            TagVar Var;

            // --- automatic resolution ---
            const VarInfo& info = registry.get(name);
            Var.name    = name;
            Var.index   = info.index;
            Var.derived = info.is_derived;

            // --- read refinement options ---
            ParmParse ppr("run.refine." + name);

            ppr.query("max_level", Var.max_level);
            ppr.queryarr("grad",          Var.grad);
            ppr.queryarr("value_less",    Var.value_less);
            ppr.queryarr("value_greater", Var.value_greater);
            ppr.queryarr("value_range",   Var.value_range);

            Var.do_grad          = !Var.grad.empty();
            Var.do_value_less    = !Var.value_less.empty();
            Var.do_value_greater = !Var.value_greater.empty();
            Var.do_value_range   = !Var.value_range.empty();

            // ErrorEst reads one threshold for each level it tags (two for value_range)
            const int n_tag_lev = std::max(0, std::min(Var.max_level, max_level));
            auto check_size = [&](amrex::Vector<amrex::Real> const& vals,
                                  int per_lev, std::string const& opt) {
                if (!vals.empty() && vals.size() < per_lev*n_tag_lev) {
                    amrex::Abort("run.refine." + name + "." + opt + " has " + std::to_string(vals.size())
                                 + " values, but " + std::to_string(per_lev*n_tag_lev)
                                 + " are needed to tag levels 0 to " + std::to_string(n_tag_lev-1)
                                 + ". Give more values or lower run.refine." + name + ".max_level.");
                }
            };
            check_size(Var.grad,          1, "grad");
            check_size(Var.value_less,    1, "value_less");
            check_size(Var.value_greater, 1, "value_greater");
            check_size(Var.value_range,   2, "value_range");

            tag_vars[i] = std::move(Var);
        }

       
    }

    ReadParameters_Derived();

#ifdef DYNAMIC_INIT
    h_parm->DynamicInit();
#endif
    h_parm->Initialize();

#if (AMREX_USE_GPU == true)
    amrex::Gpu::htod_memcpy(d_parm, h_parm, sizeof(Parm));
#else
    d_parm=h_parm;
#endif
}


void
COMPAS::AverageDownTo (Vector<MultiFab>& mf, int crse_lev)
{

    BL_PROFILE("AverageDownTo()");

    amrex::average_down(mf[crse_lev+1], mf[crse_lev],
                    geom[crse_lev+1], geom[crse_lev],
                    0, mf[crse_lev].nComp(), refRatio(crse_lev));    
}


// set covered coarse cells to be the average of overlying fine cells
void
COMPAS::AverageDown (Vector<MultiFab>& mf)
{

    BL_PROFILE("AverageDown()");

    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo (mf, lev);
    }
}

void
COMPAS::AverageDownFacesTo (Vector< Array<MultiFab,AMREX_SPACEDIM> >& mf, int crse_lev)
{

    BL_PROFILE("AverageDownFacesTo()");

    amrex::average_down_faces(amrex::GetArrOfConstPtrs(mf[crse_lev+1]),
                              amrex::GetArrOfPtrs     (mf[crse_lev  ]),
                              refRatio(crse_lev), Geom(crse_lev));
}

void
COMPAS::AverageDownFaces (Vector< Array<MultiFab,AMREX_SPACEDIM> >& mf)
{

    BL_PROFILE("AverageDownFaces()");
    
    for (int lev = finest_level-1; lev >= 0; lev--)
    {
        AverageDownFacesTo (mf, lev);
    }
}

// compute a new multifab by copying in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
COMPAS::FillPatch (MultiFab& mf,
                     int lev,
                     Real const& time,
                     Vector<MultiFab>& mf_new,
                     Vector<Real>&      t_new,
                     Vector<MultiFab>& mf_old,
                     Vector<Real>&      t_old,
                     int scomp,
                     int dcomp,
                     int ncomp)
{
    Interpolater* mapper;
    if (ID_Interpolater == "pc_interp"){
        mapper = &pc_interp;
    } else if (ID_Interpolater == "lincc_interp"){
        mapper = &lincc_interp;
    } else if (ID_Interpolater == "quartic_interp"){
        mapper = &quartic_interp;
    }

    FillPatch(mf, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              scomp, dcomp, ncomp,
              bcs, scomp,
              mapper);
}


void
COMPAS::FillPatch (MultiFab& mf,
                     int lev,
                     Real const& time,
                     Vector<MultiFab>& mf_new,
                     Vector<Real>&      t_new,
                     Vector<MultiFab>& mf_old,
                     Vector<Real>&      t_old,
                     int scomp,
                     int dcomp,
                     int ncomp,
                     Vector<BCRec> const& bcs,
                     int bcscomp,
                     Interpolater* mapper)
{

    BL_PROFILE("FillPatch[FV]()");

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(smf, stime,
                lev, time,
                mf_new, t_new,
                mf_old, t_old);
/*
        if(Gpu::inLaunchRegion())
*/
        {
            GpuBndryFuncFab<AmrCoreGPUFill> gpu_bndry_func(AmrCoreGPUFill{d_parm});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > physbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchSingleLevel(mf, time,
                                        smf, stime,
                                        scomp, dcomp, ncomp,
                                        geom[lev], physbc, bcscomp);
        }
/*
        else
        {
            AmrCoreCPUFill::parm = d_parm;
            CpuBndryFuncFab bndry_func(AmrCoreCPUFill::Boundary_Fill);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time,
                                        smf, stime,
                                        scomp, dcomp, ncomp,
                                        geom[lev], physbc, bcscomp);
        }
*/
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(cmf, ctime,
                lev-1, time,
                mf_new, t_new,
                mf_old, t_old);
        GetData(fmf, ftime,
                lev, time,
                mf_new, t_new,
                mf_old, t_old);
/*
        if(Gpu::inLaunchRegion())
*/
        {
            GpuBndryFuncFab<AmrCoreGPUFill> gpu_bndry_func(AmrCoreGPUFill{d_parm});
            PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
            PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
            amrex::FillPatchTwoLevels(mf, time,
                                      cmf, ctime,
                                      fmf, ftime,
                                      scomp, dcomp, ncomp,
                                      geom[lev-1], geom[lev],
                                      cphysbc, bcscomp, fphysbc, bcscomp,
                                      refRatio(lev-1), mapper,
                                      bcs, bcscomp);
        }
/*
        else
        {
            AmrCoreCPUFill::parm = d_parm;
            CpuBndryFuncFab bndry_func(AmrCoreCPUFill::Boundary_Fill);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchTwoLevels(mf, time,
                                      cmf, ctime,
                                      fmf, ftime,
                                      scomp, dcomp, ncomp,
                                      geom[lev-1], geom[lev],
                                      cphysbc, bcscomp, fphysbc, bcscomp,
                                      refRatio(lev-1), mapper,
                                      bcs, bcscomp);
        }
*/
    }
}


// fill mf for a Runge-Kutta stage that starts at stage_time
// FillPatch picks the fine-level state by its time label, from mf_new at t_new or mf_old at t_old,
// so it cannot fill a stage whose state is not mf_new or whose time is not t_new. Here the valid
// region and the same-level ghost cells come from the stage state mf_stage, the physical boundary
// conditions are evaluated at stage_time, and the coarse-fine ghost cells are interpolated from
// the coarser level after interpolating it in time between its old and new data to stage_time
void
COMPAS::FillPatchStage (MultiFab& mf,
                        int lev,
                        Real const& stage_time,
                        MultiFab& mf_stage,
                        Vector<MultiFab>& mf_new,
                        Vector<Real>&      t_new,
                        Vector<MultiFab>& mf_old,
                        Vector<Real>&      t_old,
                        int scomp,
                        int dcomp,
                        int ncomp)
{

    BL_PROFILE("FillPatchStage[FV]()");

    Interpolater* mapper;
    if (ID_Interpolater == "pc_interp"){
        mapper = &pc_interp;
    } else if (ID_Interpolater == "lincc_interp"){
        mapper = &lincc_interp;
    } else if (ID_Interpolater == "quartic_interp"){
        mapper = &quartic_interp;
    }

    // the stage state is the only fine-level source, labeled with the stage time
    Vector<MultiFab*> fmf{&mf_stage};
    Vector<Real> ftime{stage_time};

    GpuBndryFuncFab<AmrCoreGPUFill> gpu_bndry_func(AmrCoreGPUFill{d_parm});

    if (lev == 0)
    {
        PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > physbc(geom[lev],bcs,gpu_bndry_func);
        amrex::FillPatchSingleLevel(mf, stage_time,
                                    fmf, ftime,
                                    scomp, dcomp, ncomp,
                                    geom[lev], physbc, scomp);
    }
    else
    {
        Vector<MultiFab*> cmf;
        Vector<Real> ctime;
        GetData(cmf, ctime,
                lev-1, stage_time,
                mf_new, t_new,
                mf_old, t_old);

        PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
        amrex::FillPatchTwoLevels(mf, stage_time,
                                  cmf, ctime,
                                  fmf, ftime,
                                  scomp, dcomp, ncomp,
                                  geom[lev-1], geom[lev],
                                  cphysbc, scomp, fphysbc, scomp,
                                  refRatio(lev-1), mapper,
                                  bcs, scomp);
    }
}


// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
COMPAS::FillCoarsePatch (MultiFab& mf,
                           int lev,
                           Real const& time,
                           Vector<MultiFab>& mf_new,
                           Vector<Real>&      t_new,
                           Vector<MultiFab>& mf_old,
                           Vector<Real>&      t_old,
                           int scomp,
                           int dcomp,
                           int ncomp)
{

    Interpolater* mapper;
    if (ID_Interpolater == "pc_interp"){
        mapper = &pc_interp;
    } else if (ID_Interpolater == "lincc_interp"){
        mapper = &lincc_interp;
    } else if (ID_Interpolater == "quartic_interp"){
        mapper = &quartic_interp;
    }

    FillCoarsePatch(mf, lev, time,
                    mf_new, t_new,
                    mf_old, t_old,
                    scomp, dcomp, ncomp,
                    bcs, scomp,
                    mapper);
}


void
COMPAS::FillCoarsePatch (MultiFab& mf,
                           int lev,
                           Real const& time,
                           Vector<MultiFab>& mf_new,
                           Vector<Real>&      t_new,
                           Vector<MultiFab>& mf_old,
                           Vector<Real>&      t_old,
                           int scomp,
                           int dcomp,
                           int ncomp,
                           Vector<BCRec> const& bcs,
                           int bcscomp,
                           Interpolater* mapper)
{

    BL_PROFILE("FillCoarsePatch[FV]()");

    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(cmf, ctime,
            lev-1, time,
            mf_new, t_new,
            mf_old, t_old);

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }
/*
    if(Gpu::inLaunchRegion())
*/
    {
        GpuBndryFuncFab<AmrCoreGPUFill> gpu_bndry_func(AmrCoreGPUFill{d_parm});
        PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > cphysbc(geom[lev-1],bcs,gpu_bndry_func);
        PhysBCFunct<GpuBndryFuncFab<AmrCoreGPUFill> > fphysbc(geom[lev],bcs,gpu_bndry_func);
        amrex::InterpFromCoarseLevel(mf, time,
                                     *cmf[0],
                                     scomp, dcomp, ncomp,
                                     geom[lev-1], geom[lev],
                                     cphysbc, bcscomp, fphysbc, bcscomp,
                                     refRatio(lev-1), mapper,
                                     bcs, bcscomp);
    }
/*
    else
    {
        AmrCoreCPUFill::parm = d_parm;
        CpuBndryFuncFab bndry_func(AmrCoreCPUFill::Boundary_Fill);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);
        amrex::InterpFromCoarseLevel(mf, time,
                                     *cmf[0],
                                     scomp, dcomp, ncomp,
                                     geom[lev-1], geom[lev],
                                     cphysbc, bcscomp, fphysbc, bcscomp,
                                     refRatio(lev-1), mapper,
                                     bcs, bcscomp);
    }
*/
}


// utility to copy in data from dof_old and/or dof_new into another multifab
void
COMPAS::GetData (Vector<MultiFab*>& data,
                   Vector<Real>& datatime,
                   int lev,
                   Real const& time,
                   Vector<MultiFab>& mf_new,
                   Vector<Real>&      t_new,
                   Vector<MultiFab>& mf_old,
                   Vector<Real>&      t_old)
{
    /* CHANGE MF HERE */
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.0e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(&mf_new[lev]);
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(&mf_old[lev]);
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(&mf_old[lev]);
        data.push_back(&mf_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

void
COMPAS::RefluxLev (int lev)
{
#if (NONCONSERVATIVE == true)
    amrex::Abort("Need to define a nonconservative reflux procedure!\n");
#endif
    // update lev based on coarse-fine flux mismatch
    flux_reg[lev+1]->Reflux(dof_new[lev], 1.0, 0, 0, NSTATE, geom[lev]);
}

void
COMPAS::AdvanceAllLevels (Real time, Real dt_lev, int ncycle)
{
    for (int lev = 0; lev <= finest_level; lev++){
        t_old[lev] = t_new[lev];
        t_new[lev] += dt_lev;
    }
    
    for (int lev = 0; lev <= finest_level; lev++){
        if (ID_TimeIntegrator == "ForwardEuler"){
            ForwardEuler (dof_new,dof_old,lev,time,dt_lev,ncycle);
        }  else if (ID_TimeIntegrator == "TVD-RK2"){
            SecondOrderSSPRK (dof_new,dof_old,lev,time,dt_lev,ncycle);
        } else if (ID_TimeIntegrator == "TVD-RK3"){
            ThirdOrderSSPRK (dof_new,dof_old,lev,time,dt_lev,ncycle);
        } else if (ID_TimeIntegrator == "RK4"){
            FourthOrderRK (dof_new,dof_old,lev,time,dt_lev,ncycle);
        }
    }

    if (do_reflux){
        for (int lev = 0; lev < finest_level; lev++){
            RefluxLev(lev);
        }
    }   

}

void 
COMPAS::AdvanceAtLevel (Vector<MultiFab>& mf_new, 
                         Vector<MultiFab>& mf_old, 
                         int lev, 
                         Real time, 
                         Real dt_lev, 
                         int ncycle)
{

    if (ID_TimeIntegrator == "ForwardEuler"){
        ForwardEuler (mf_new,mf_old,lev,time,dt_lev,ncycle);
    }  else if (ID_TimeIntegrator == "TVD-RK2"){
        SecondOrderSSPRK (mf_new,mf_old,lev,time,dt_lev,ncycle);
    } else if (ID_TimeIntegrator == "TVD-RK3"){
        ThirdOrderSSPRK (mf_new,mf_old,lev,time,dt_lev,ncycle);
    } else if (ID_TimeIntegrator == "RK4"){
        FourthOrderRK (mf_new,mf_old,lev,time,dt_lev,ncycle);
    }

}

// Advance a level by dt
// (includes a recursive call for finer levels)
void
COMPAS::timeStepWithSubcycling (Vector<MultiFab>& mf_new,
                                  Vector<MultiFab>& mf_old,
                                  int lev, 
                                  amrex::Real time, 
                                  int iteration,
                                  int stage)
{

    BL_PROFILE("timeStepWithSubcycling()");

    if (stage == 0){
        if (regrid_int > 0)  // We may need to regrid
        {

            // help keep track of whether a level was already regridded
            // from a coarser level call to regrid
            static Vector<int> last_regrid_step(max_level+1, 0);

            // regrid changes level "lev+1" so we don't regrid on max_level
            // also make sure we don't regrid fine levels again if
            // it was taken care of during a coarser regrid
            if (lev < max_level && istep[lev] > last_regrid_step[lev])
            {
                if (istep[lev] % regrid_int == 0)
                {
                    // regrid could add newly refine levels (if finest_level < max_level)
                    // so we save the previous finest level index
                    int old_finest = finest_level;
                    regrid(lev, time);

                    // mark that we have regridded this level already
                    for (int k = lev; k <= finest_level; ++k) {
                        last_regrid_step[k] = istep[k];
                    }

                    // if there are newly created levels, set the time step
                    for (int k = old_finest+1; k <= finest_level; ++k) {
#if (ADVECTION == true && DIFFUSION == true)
                        dt[k] = dt[k-1] / (MaxRefRatio(k-1)*MaxRefRatio(k-1));
#else
                        dt[k] = dt[k-1] / MaxRefRatio(k-1);
#endif
                    }
                }
            }
        }
    }

    if (stage == 0){
        if (Verbose()) {
            amrex::Print() << std::string(4*lev, ' ') << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << std::string(4*(finest_level - lev) + 2, ' ') << "Advance with time = " << t_new[lev]
                           << " dt = " << dt[lev];
        }
    }

    t_old[lev] = t_new[lev];
    t_new[lev] += dt[lev];

    // Advance a single level for a single time step, and update flux registers
    AdvanceAtLevel(mf_new, mf_old, lev, time, dt[lev], nsubsteps[lev]);

    if (stage == 0){
        ++istep[lev];

        if (Verbose())
        {
            // amrex::Print() << std::string(4*lev, ' ') << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << " ==> Advanced " << CountCells(lev) << " cells" << std::endl;
        }
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStepWithSubcycling(mf_new, mf_old, lev+1, time+(i-1)*dt[lev+1], i, stage);
        }

        if (do_reflux)
        {
            RefluxLev(lev);
        }

        AverageDownTo(mf_new,lev); // average lev+1 down to lev
    }

}

// Advance all the levels with the same dt
void
COMPAS::timeStepNoSubcycling (Real time, int iteration)
{
    if (max_level > 0 && regrid_int > 0)  // We may need to regrid
    {
        if (istep[0] % regrid_int == 0)
        {
            regrid(0, time);
        }
    }

    if (Verbose()) {
        for (int lev = 0; lev <= finest_level; lev++)
        {
           amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
           amrex::Print() << "ADVANCE with time = " << t_new[lev]
                          << " dt = " << dt[lev] << std::endl;
        }
    }

    AdvanceAllLevels (time, dt[0], iteration);

    // Make sure the coarser levels are consistent with the finer levels
    AverageDown (dof_new);

    for (int lev = 0; lev <= finest_level; lev++)
        ++istep[lev];

    if (Verbose())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }
    }

}


void
COMPAS::FixedDt ()
{

    amrex::Real dt_0 = dt_fixed/n_cells[0];

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;

    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;

    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }

}


void
COMPAS::OutputError ()
{

    amrex::Vector<amrex::Real> Linf(NSTATE);
    amrex::Vector<amrex::Real> L2(NSTATE);

    EvaluateError(L2, Linf);

    for (int v = 0; v < NSTATE; v++){
        ParallelDescriptor::ReduceRealSum(&L2[v],  1); // Get sum of L2 errors across processors
        ParallelDescriptor::ReduceRealMax(&Linf[v],1); // Get max Linf error across processors
    }

    for (int v = 0; v < NSTATE; v++){
        L2[v] = std::sqrt(L2[v]);
    }

    int myProc = ParallelDescriptor::MyProc();

    if (myProc == 0){
        /* Write the results and the simulation parameters to the specified file */
        std::ofstream ofs;
        ofs.open(conv_output_file + ".0", std::ofstream::out | std::ofstream::trunc);
        ofs.close();
        amrex::PrintToFile fp(conv_output_file);
        fp.SetPrecision(16);
        // fp << "ProblemName PhysicsName SpatialScheme TemporalScheme Nx Ny dt Runtime L2Error LinfError\n";
        fp << conv_prob_name  << " " 
           << n_cells[0] << " " 
           << n_cells[1] << " " 
           << dt_fixed << " " 
           << Evolve_stop_time << "\n";

        for (int v = 0; v < NSTATE; v++){
            fp << vnamesC[v] << " "
               << L2[v] << " " 
               << Linf[v] << "\n";
        } 
    }
           

}



// a wrapper for EstTimeStep
void
COMPAS::ComputeDt ()
{

    BL_PROFILE("ComputeDt()");

    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EstTimeStep(dt_tmp[lev], lev, t_new[lev]);
        if (t_new[lev] == 0.0){
            dt[lev] = dt_tmp[lev];
        }
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;

    for (int lev = 0; lev <= finest_level; ++lev) {
        if (timestep_change_limiter && dt[lev] > 0.0){
            dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        }
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = (1.0e-3)*dt_0;

    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;

    if (do_subcycle){
        for (int lev = 1; lev <= finest_level; ++lev) {
            dt[lev] = dt[lev-1] / nsubsteps[lev];
        }
    }
    else{
        amrex::Real dt_min = dt_0;
        for (int lev = 0; lev <= finest_level; ++lev) {
            dt_min = std::min(dt_min, dt_tmp[lev]);
        }
        for (int lev = 0; lev <= finest_level; ++lev) {
            dt[lev] = dt_min;
        }
    }

}

// compute dt from CFL considerations
void
COMPAS::EstTimeStep (Real& dt_est, int lev, Real const& time)
{
    BL_PROFILE("COMPAS::EstTimeStep()");

    dt_est = std::numeric_limits<Real>::max();

    const Real* dx  =  geom[lev].CellSize();

    if (time == 0.0) {
       CalculateInitialcmax(lev);
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        
#if (ADVECTION == true)
        // Real est = facevel[lev][idim].norminf(0,0,true);
        Real est_adv = c_max[lev].norminf(0,0,true);
        dt_est = amrex::min(dt_est, cfl*dx[idim]/est_adv);
#endif

#if (DIFFUSION==true)
        Real est_diff = c_max[lev].norminf(1,0,true);
        dt_est = amrex::min(dt_est, vnn*dx[idim]*dx[idim]/est_diff);
#endif
    }
    // dt_est *= cfl;
    /* TODO: Test here if the above line matters or we can pre-mulitply cfl and vnn */
}


void 
COMPAS::variableSetUp()
{

    lo_bc.resize(3); hi_bc.resize(3);

    h_parm = new Parm{}; 
    d_parm = (Parm*)The_Arena()->alloc(sizeof(Parm));

}


void 
COMPAS::variableCleanUp()
{
    // Doesn't seem to need this, keep it around anyway for now
    delete h_parm;
    The_Arena()->free(d_parm);
}


// get plotfile name
std::string
COMPAS::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}



void
COMPAS::ResetOutputMF ()
{
    BL_PROFILE("ResetOutputMF()");
    
    if (num_output_vars != NSTATE){
        for (int lev = 0; lev <= finest_level; lev++){
            dof_old[lev].define(grids[lev], dmap[lev], NSTATE, 0);
        }
    }


    
}

void
COMPAS::CalculateOutputMF ()
{

    BL_PROFILE("CalculateOutputMF()");

    if (output_state_vars){
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::Copy(dof_old[lev],dof_new[lev],0,0,NSTATE,0);
        }
    } else {
        if (num_output_vars != NSTATE){
            for (int lev = 0; lev <= finest_level; ++lev) {
                dof_old[lev].define(grids[lev], dmap[lev], num_output_vars, 0);    
            }
        }

        // The derived variables from the velocity gradient need the neighbors of each cell
        bool need_ghosts = false;
        for (const std::string& name : output_vars) {
            need_ghosts = need_ghosts || DerivedVarNeedsNeighbors(name);
        }

        for (int lev = 0; lev <= finest_level; ++lev) {

            MultiFab& state_old = dof_old[lev];
            MultiFab Sborder;
            if (need_ghosts) {
                // dof_old holds the output now, so dof_new is passed for both times
                Sborder.define(grids[lev], dmap[lev], NSTATE, 1);
                FillPatch(Sborder, lev, t_new[lev],
                          dof_new, t_new,
                          dof_new, t_new,
                          0, 0, NSTATE);
            }
            MultiFab& state_new = need_ghosts ? Sborder : dof_new[lev];

            // MultiFab OutputMF(grids[lev], dmap[lev], num_output_vars, 0);

            const auto problo = Geom(lev).ProbLoArray();
            const auto dx_lev = Geom(lev).CellSizeArray();
            // Cell sizes in all three directions, 1 in the directions that do not exist
            amrex::GpuArray<Real, 3> dx = {1.0, 1.0, 1.0};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) dx[d] = dx_lev[d];

            Parm const* lparm = d_parm;
            amrex::Array<int, 2*NSTATE> const out_vars_t = out_vars;
            const int nvars = num_output_vars;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(state_new,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Array4<Real      > const output_vars = state_old.array(mfi);
                Array4<Real const> const state_vars  = state_new.array(mfi);
                const Box& box = mfi.tilebox();

                amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Get derived variables here and fill them in dof_old
                    compute_derived_var (i,j,k,out_vars_t,nvars, dx[0], dx[1], dx[2], output_vars, state_vars, *lparm);
                });

            }
            // std::swap(OutputMF, dof_old[lev]);
        }
    }
}





// put together an array of multifabs for writing
Vector<const MultiFab*>
COMPAS::PlotFileMF () const
{

    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&dof_old[i]);    
        
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
COMPAS::PlotFileVarNames () const
{
    if (output_state_vars){
        return vnamesC;
    } else{
        return output_vars;
    }
}

// write plotfile to disk
void
COMPAS::WritePlotFile () const
{

    BL_PROFILE("WritePlotFile()");

    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), t_new[0], istep, refRatio());


    const char *movie_file_c = movie_filename.c_str();
    std::string writename = "plt_"+case_name;
    const char *wname = writename.c_str();
    std::fstream movie_file;
    if (istep[0] == 0){
        movie_file.open(movie_file_c, std::ofstream::out | std::ofstream::trunc);
    } else {
        movie_file.open(movie_file_c, std::ios::app);
    }
    if (amrex::ParallelDescriptor::MyProc() == 0){
        movie_file << amrex::Concatenate(wname, istep[0], 5) + "/Header" << "\n";
    }
    movie_file.close();

}

void
COMPAS::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for COMPAS\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   // The header version without per-box min/max is used because subnormal min/max values
   // (e.g. a volume fraction near zero) cannot be read back on some platforms at restart.
   VisMF::Header::Version header_version = VisMF::GetHeaderVersion();
   VisMF::SetHeaderVersion(VisMF::Header::NoFabHeader_v1);
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(dof_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", ""));
   }
   VisMF::SetHeaderVersion(header_version);

}

namespace {
// utility to skip to next line in Header
void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}


void
COMPAS::Coarsen (int ref_ratio)
{
    amrex::IntVect tempRefRatio(ref_ratio);
    for (int lev = 0; lev <= finest_level; lev++){
        geom[lev].coarsen(tempRefRatio);
        grids[lev].coarsen(ref_ratio);
    }
}

void
COMPAS::Coarsen (int ref_ratio, int lev)
{
    amrex::IntVect tempRefRatio(ref_ratio);
    geom[lev].coarsen(tempRefRatio);
    grids[lev].coarsen(ref_ratio);
}

void
COMPAS::Refine (int ref_ratio)
{
    amrex::IntVect tempRefRatio(ref_ratio);
    for (int lev = 0; lev <= finest_level; lev++){
        geom[lev].refine(tempRefRatio);
        grids[lev].refine(ref_ratio);    
    }
    
}

void
COMPAS::Refine (int ref_ratio, int lev)
{
    amrex::IntVect tempRefRatio(ref_ratio);
    geom[lev].refine(tempRefRatio);
    grids[lev].refine(ref_ratio);
}


void
COMPAS::ReadCheckpointFile ()
{

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = NSTATE;
        int ncomp_nc = NC_TERMS; // Change if needed
        int nghost = NGROW;
        dof_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        dof_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

        if (ID_IB){
            ib_levelset[lev].define(grids[lev], dmap[lev], 1, nghost);
        }

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, FR_COMP));
#if (NONCONSERVATIVE == true)
            flux_reg_nc[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp_nc));
#endif
        }

        int n_cmax_comp = (ADVECTION == true) + (DIFFUSION == true);
        c_max[lev].define(ba, dm, n_cmax_comp, 1);

        MakeNewLevelFromScratch_Derived(lev, t_new[0], ba, dm);
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(dof_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", ""));
    }

}

void
COMPAS::SaveLine (amrex::Vector<amrex::MultiFab> & mf, 
                    int dir, 
                    amrex::Array<amrex::Real,2> dir_coord, 
                    amrex::Vector<std::string> vnames,
                    std::string fname)
{

    amrex::Vector<amrex::MultiFab> SliceData(max_level + 1);

    int dir1, dir2;
    if (dir == 0){
        dir1 = 1; dir2 = 2;
    } else if (dir == 1){
        dir1 = 0; dir2 = 2;
    } else if (dir == 2){
        dir1 = 0; dir2 = 1;
    }

    for (int lev = 0; lev < finest_level + 1; lev++){

        const Geometry& geom_full = geom[lev];
        RealBox real_slice = geom_full.ProbDomain();
        real_slice.setLo(dir1, dir_coord[0]);
        real_slice.setHi(dir1, dir_coord[0] + geom_full.CellSize(dir1));
#if AMREX_SPACEDIM > 2
        real_slice.setLo(dir2, dir_coord[1]);
        real_slice.setHi(dir2, dir_coord[1] + geom_full.CellSize(dir2));
#endif

        IntVect slice_lo, slice_hi;
        AMREX_D_TERM(slice_lo[0]=static_cast<int>(std::floor((real_slice.lo(0) - geom_full.ProbLo(0))/geom_full.CellSize(0)));,
                   slice_lo[1]=static_cast<int>(std::floor((real_slice.lo(1) - geom_full.ProbLo(1))/geom_full.CellSize(1)));,
                   slice_lo[2]=static_cast<int>(std::floor((real_slice.lo(2) - geom_full.ProbLo(2))/geom_full.CellSize(2))););
        AMREX_D_TERM(slice_hi[0]=static_cast<int>(std::floor((real_slice.hi(0) - geom_full.ProbLo(0))/geom_full.CellSize(0)));,
                   slice_hi[1]=static_cast<int>(std::floor((real_slice.hi(1) - geom_full.ProbLo(1))/geom_full.CellSize(1)));,
                   slice_hi[2]=static_cast<int>(std::floor((real_slice.hi(2) - geom_full.ProbLo(2))/geom_full.CellSize(2))););

        slice_hi[dir1] = slice_lo[dir1]+1;
#if AMREX_SPACEDIM > 2 
        slice_hi[dir2] = slice_lo[dir2]+1;
#endif

        Box slice_box(slice_lo, slice_hi);
        BoxArray slice_ba(slice_box);
        slice_ba = intersect(slice_ba, mf[lev].boxArray());
        DistributionMapping slice_dm(slice_ba);
        SliceData[lev].define(slice_ba, slice_dm, mf[0].nComp(), 0);
        SliceData[lev].ParallelCopy(mf[lev]);
    }

    amrex::WriteMultiLevelPlotfile(amrex::Concatenate(fname, istep[0], 5),
                                   finest_level+1, 
                                   amrex::GetVecOfConstPtrs(SliceData), 
                                   vnames,
                                   Geom(), t_new[0], istep, refRatio());


}


void
COMPAS::SavePlane (amrex::Vector<amrex::MultiFab> & mf, 
                     int dir, 
                     amrex::Real dir_coord, 
                     amrex::Vector<std::string> vnames,
                     std::string fname)
{
    amrex::Vector<amrex::MultiFab> SliceData(max_level + 1);

    for (int lev = 0; lev < finest_level + 1; lev++){

        const Geometry& geom_full = geom[lev];
        RealBox real_slice = geom_full.ProbDomain();
        real_slice.setLo(dir, dir_coord);
        real_slice.setHi(dir, dir_coord + geom_full.CellSize(dir));

        IntVect slice_lo, slice_hi;
        AMREX_D_TERM(slice_lo[0]=static_cast<int>(std::floor((real_slice.lo(0) - geom_full.ProbLo(0))/geom_full.CellSize(0)));,
                   slice_lo[1]=static_cast<int>(std::floor((real_slice.lo(1) - geom_full.ProbLo(1))/geom_full.CellSize(1)));,
                   slice_lo[2]=static_cast<int>(std::floor((real_slice.lo(2) - geom_full.ProbLo(2))/geom_full.CellSize(2))););
        AMREX_D_TERM(slice_hi[0]=static_cast<int>(std::floor((real_slice.hi(0) - geom_full.ProbLo(0))/geom_full.CellSize(0)));,
                   slice_hi[1]=static_cast<int>(std::floor((real_slice.hi(1) - geom_full.ProbLo(1))/geom_full.CellSize(1)));,
                   slice_hi[2]=static_cast<int>(std::floor((real_slice.hi(2) - geom_full.ProbLo(2))/geom_full.CellSize(2))););

        slice_hi[dir] = slice_lo[dir]+1;
        Box slice_box(slice_lo, slice_hi);
        BoxArray slice_ba(slice_box);
        slice_ba = intersect(slice_ba, mf[lev].boxArray());
        DistributionMapping slice_dm(slice_ba);
        SliceData[lev].define(slice_ba, slice_dm, mf[lev].nComp(), 0);
        SliceData[lev].ParallelCopy(mf[lev]);
    }

    amrex::WriteMultiLevelPlotfile(amrex::Concatenate(fname, istep[0], 5),
                                   finest_level+1, 
                                   amrex::GetVecOfConstPtrs(SliceData), 
                                   vnames,
                                   Geom(), t_new[0], istep, refRatio());

}

void
COMPAS::EvaluateError (amrex::Vector<amrex::Real> & L2, 
                         amrex::Vector<amrex::Real> & Linf)
{
#if (CONVERGENCE==true)
    /* SINGLE LEVEL FOR NOW */
    int lev = 0; 

    const auto          prob_lo        = Geom(lev).ProbLoArray();
    const auto          dx             = Geom(lev).CellSizeArray();
    const auto          dx_FinestLevel = dxFinest();
    const amrex::Real   t_lev          = t_new[0];

    Parm const* lparm = d_parm;

    
    for (int v = 0; v < NSTATE; v++){
        L2[v] = 0.0;
        Linf[v] = 0.0;
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    FArrayBox tmpfab;

    for (MFIter mfi(dof_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();
        Array4<Real> statein   = dof_new[lev].array(mfi);

        tmpfab.resize(bx,2*NSTATE);
        tmpfab.setVal<RunOn::Device>(Real(0.0));
        Elixir tmpeli = tmpfab.elixir();
        Array4<Real> err_arr = tmpfab.array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            using ErrQuad = Quadrature::GaussLegendre4;
            constexpr int nq = ErrQuad::npts;
            constexpr amrex::Real ref_vol_inv = 1.0_rt / (AMREX_D_TERM(2.0_rt,*2.0_rt,*2.0_rt));
            const amrex::Real cell_vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
            amrex::Array<amrex::Real, NSTATE> Conservative_Exact_Avg;

            for (int i_State = 0; i_State < NSTATE; i_State++) {
                Conservative_Exact_Avg[i_State] = 0.0_rt;
            }

            for (int qx = 0; qx < nq; ++qx)
            {
                const amrex::Real x = prob_lo[0] +
                    (static_cast<amrex::Real>(i) + 0.5_rt + 0.5_rt*ErrQuad::node(qx)) * dx[0];
                const amrex::Real wx = ErrQuad::weight(qx);

#if (AMREX_SPACEDIM > 1)
                for (int qy = 0; qy < nq; ++qy)
                {
                    const amrex::Real y = prob_lo[1] +
                        (static_cast<amrex::Real>(j) + 0.5_rt + 0.5_rt*ErrQuad::node(qy)) * dx[1];
                    const amrex::Real wy = ErrQuad::weight(qy);

#if (AMREX_SPACEDIM > 2)
                    for (int qz = 0; qz < nq; ++qz)
                    {
                        const amrex::Real z = prob_lo[2] +
                            (static_cast<amrex::Real>(k) + 0.5_rt + 0.5_rt*ErrQuad::node(qz)) * dx[2];
                        const amrex::Real wz = ErrQuad::weight(qz);
                        const amrex::Real wq = wx * wy * wz;
#else
                    {
                        const amrex::Real z = 0.0_rt;
                        const amrex::Real wq = wx * wy;
#endif
                        amrex::Array<amrex::Real, NSTATE> Conservative_Exact_Point;
                        exact_solution(Conservative_Exact_Point, x, y, z, t_lev,
                            AMREX_D_DECL(dx[0],dx[1],dx[2]),
                            AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                            *lparm);

                        for (int i_State = 0; i_State < NSTATE; i_State++)
                        {
                            Conservative_Exact_Avg[i_State] += wq * Conservative_Exact_Point[i_State];
                        }
                    }
#if (AMREX_SPACEDIM > 2)
                    }
#endif
                }
#else
                {
                    const amrex::Real y = 0.0_rt;
                    const amrex::Real z = 0.0_rt;
                    const amrex::Real wq = wx;

                    amrex::Array<amrex::Real, NSTATE> Conservative_Exact_Point;
                    exact_solution(Conservative_Exact_Point, x, y, z, t_lev,
                        AMREX_D_DECL(dx[0],dx[1],dx[2]),
                        AMREX_D_DECL(dx_FinestLevel[0],dx_FinestLevel[1],dx_FinestLevel[2]),
                        *lparm);

                    for (int i_State = 0; i_State < NSTATE; i_State++)
                    {
                        Conservative_Exact_Avg[i_State] += wq * Conservative_Exact_Point[i_State];
                    }
                }
#endif
            }

            for (int i_State = 0; i_State < NSTATE; i_State++) {
                Conservative_Exact_Avg[i_State] *= ref_vol_inv;

                const amrex::Real diff =
                    statein(i,j,k,i_State) - Conservative_Exact_Avg[i_State];
                err_arr(i,j,k,i_State) = diff * diff * cell_vol;
                err_arr(i,j,k,i_State + NSTATE) = std::abs(diff);
            }

        });

        for (int i_State = 0; i_State < NSTATE; i_State++){
            L2[i_State]  += tmpfab.sum<RunOn::Device>(i_State);
            Linf[i_State] = std::max({Linf[i_State], tmpfab.max<RunOn::Device>(i_State + NSTATE)});
        }
    }
#endif
}















