// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
// Portions derived from AMReX, BSD-3-Clause. See THIRD_PARTY_NOTICES.md.

#include <COMPAS.H>
#include <Kernels.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Gpu.H>
#include <LinearAlgebra.H>
#include <Tools.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <Class_Compressible_PhaseField.H>

#include <Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==FIVEEQS_NPHASE || PHYSICS==SIXEQS || PHYSICS==SIXEQS_IE_NPHASE)
#include <IncludePhysics.H>
#include <Parm.H>
#include <Operators_LinearSystem.H>
#include <Physics_PhaseField.H>
#include <Operators_PhaseField.H>

#if (PHYSICS==FIVEEQS || PHYSICS==SIXEQS)
#include <Physics_Compressible_TwoPhase.H>
#include <Physics_Compressible_TwoPhase_PhaseField.H>
#endif

#if (PHYSICS==FIVEEQS_NPHASE || PHYSICS==SIXEQS_IE_NPHASE)
#include <Physics_Compressible_NPhase.H>
#include <Physics_Compressible_NPhase_PhaseField.H>
#endif

#endif

#if (PHYSICS==FIVEEQS)
#include <Physics_Compressible5Eq.H>
#include <Operators_Compressible5Eq.H>
#include <Physics_Compressible5Eq_PhaseField.H>
#include <Operators_Compressible5Eq_PhaseField.H>
#endif

#if (PHYSICS==FIVEEQS_NPHASE)
#include <Physics_Compressible5Eq_NPhase.H>
#include <Operators_Compressible5Eq_NPhase.H>
#include <Physics_Compressible5Eq_NPhase_PhaseField.H>
#include <Operators_Compressible5Eq_NPhase_PhaseField.H>
#endif

#if (PHYSICS==SIXEQS)
#include <Physics_Compressible6Eq.H>
#include <Operators_Compressible6Eq.H>
#include <Physics_Compressible6Eq_PhaseField.H>
#include <Operators_Compressible6Eq_PhaseField.H>
#endif

#if (PHYSICS==SIXEQS_IE_NPHASE)
#include <Physics_Compressible6Eq_IE_NPhase.H>
#include <Operators_Compressible6Eq_IE_NPhase.H>
#include <Physics_Compressible6Eq_IE_NPhase_PhaseField.H>
#include <Operators_Compressible6Eq_IE_NPhase_PhaseField.H>
#endif



Compressible_PhaseField::~Compressible_PhaseField ()
{
}

void 
Compressible_PhaseField::Physics_Derived_Constructor ()
{

    ReadParameters_Derived();
    SetPhysicsBC();

    if (h_parm->PhaseField_Parm.phase_field){
        int nlevs_max = max_level + 1;
        PF_c_max.resize(nlevs_max);
    }
        
}


void 
Compressible_PhaseField::ReadParameters_Derived ()
{
    {
        std::string RiemannSolverString = "HLLC";
        ParmParse pp("Physics");
        pp.query("RiemannSolver",RiemannSolverString);
        if (RiemannSolverString == "LLF"){
            h_parm->Physics_Parm.RiemannSolver = Riemann::LLF;
            h_parm->Physics_Parm.ID_Riemann = 0;
        }
        else if (RiemannSolverString == "HLL"){
            h_parm->Physics_Parm.RiemannSolver = Riemann::HLL;
            h_parm->Physics_Parm.ID_Riemann = 1;
        }
        else if (RiemannSolverString == "HLLC"){
            h_parm->Physics_Parm.RiemannSolver = Riemann::HLLC;
            h_parm->Physics_Parm.ID_Riemann = 2;
        }
        else {
            amrex::Abort("Physics.RiemannSolver = '" + RiemannSolverString + "' is not recognized. "
                         "Valid values are LLF, HLL, HLLC.");
        }
    }

    NGROW = h_parm->FiniteVolume_Parm.Num_Grow;
}

void 
Compressible_PhaseField::MakeNewLevelFromCoarse_Derived (int lev, 
                                                         Real time, 
                                                         const BoxArray& ba,
                                                         const DistributionMapping& dm)
{
    if (h_parm->PhaseField_Parm.phase_field){
        PF_c_max[lev].define(ba, dm, 1, 1);    
    }
}

void 
Compressible_PhaseField::MakeNewLevelFromScratch_Derived (int lev, 
                                                          Real time, 
                                                          const BoxArray& ba,
                                                          const DistributionMapping& dm)
{
    if (h_parm->PhaseField_Parm.phase_field){
        PF_c_max[lev].define(ba, dm, 1, 1);    
    }
}

void 
Compressible_PhaseField::RemakeLevel_Derived (int lev, 
                                              Real time, 
                                              const BoxArray& ba,
                                              const DistributionMapping& dm)
{
    if (h_parm->PhaseField_Parm.phase_field){
        PF_c_max[lev] = MultiFab(ba, dm, 1, 1);
    }
}


void 
Compressible_PhaseField::RemakeAllLevelsFromData_Derived (int lev,
                                                          const BoxArray& ba,
                                                          const DistributionMapping& dm)
{
    if (h_parm->PhaseField_Parm.phase_field){
        PF_c_max[lev] = MultiFab(ba, dm, 1, 1);
    }
}


void 
Compressible_PhaseField::ClearLevel_Derived (int lev)
{
    if (h_parm->PhaseField_Parm.phase_field){
        PF_c_max[lev].clear();
    }
}

// advance solution to final time
void
Compressible_PhaseField::Evolve ()
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

    //Store temporary dof_new and t_new
    Vector<MultiFab> dof_temp;
    Vector<Real>       t_temp;
    Vector<Real>      dt_temp;
    Vector<int>    istep_temp;
    int     finest_level_temp;

#ifdef USER_OUTPUT_FUNC 
    UserOutputFunction(*this);
#endif


    /* Loop over time steps, stopping when the final time is reached */
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {


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

        amrex::Print() << "\n==== Coarse step " << step+1 << " ====" << std::endl;


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

        //Copy dof_new and t_new to dof_temp and t_temp
        if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
            dof_temp.resize(finest_level+1);
              t_temp.resize(finest_level+1);
             dt_temp.resize(finest_level+1);
          istep_temp.resize(finest_level+1);
   finest_level_temp = finest_level;
            for (int lev = 0; lev <= finest_level; lev++){
                const BoxArray& ba = dof_new[lev].boxArray();
                const DistributionMapping& dm = dof_new[lev].DistributionMap();
                int Num_Comp = dof_new[lev].nComp();
                int Num_Grow = dof_new[lev].nGrow();
                dof_temp[lev].define(ba, dm, Num_Comp, Num_Grow);
                MultiFab::Copy(dof_temp[lev], dof_new[lev], 0, 0, Num_Comp, Num_Grow);
                t_temp[lev] = t_new[lev];
                dt_temp[lev] = dt[lev];
                istep_temp[lev] = istep[lev];
            }
        }

        h_parm->PhaseField_Parm.Mobility = 1.0;
        Parm* lparm = d_parm;
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
    {
        lparm->PhaseField_Parm.Mobility = 1.0;
    });
        // device_modify([=] AMREX_GPU_DEVICE (){d_parm->PhaseField_Parm.Mobility = 1.0;});

        ExplicitTimeStep(cur_time);
        if (h_parm->FiniteVolume_Parm.ID_Bound == 1 && h_parm->FiniteVolume_Parm.ID_Return == 1){
            RemakeAllLevelsFromData(dof_temp, t_temp, dt_temp, istep_temp, finest_level_temp);
            amrex::Real min_EnergyBound;
            MinEnergyBound(min_EnergyBound,
                           dof_new,
                           *h_parm);
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        dof_new,
                        *h_parm);
            //Print()<<"\n"<<"min_EnergyBound before entering the time step is "<<min_EnergyBound<<" >= 0, and min_Pressure is "<<min_Pressure<<" >= "<<d_parm->p_min<<"\n";
            Print()<<"\n"<<"min_EnergyBound before entering the time step is "<<min_EnergyBound<<" >= 0, and min_Pressure is "<<min_Pressure<<"\n";
            //Recompute a time step

            // The flux kernels read ID_Lambda from the device copy
            h_parm->FiniteVolume_Parm.ID_Lambda = 0; //0: General; 1: Special (like all gases or single-phase flows)
            // device_modify([=] AMREX_GPU_DEVICE (){d_parm->FiniteVolume_Parm.ID_Lambda = 0;});
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
            {
                lparm->FiniteVolume_Parm.ID_Lambda = 0;
            });

#if (FIXED_DT == true)
            amrex::Abort("The time step is too large to preserve bound");
#else
            ComputeDt();
#endif
                
            h_parm->PhaseField_Parm.Mobility = 0.0;
            // device_modify([=] AMREX_GPU_DEVICE (){d_parm->PhaseField_Parm.Mobility = 0.0;});

            for (int lev = 0; lev <= finest_level; lev++){
                h_parm->PhaseField_Parm.Mobility = std::max( h_parm->PhaseField_Parm.Mobility, c_max[lev].norminf(0,0,true) );
                // device_modify([=] AMREX_GPU_DEVICE (){d_parm->PhaseField_Parm.Mobility = h_parm->PhaseField_Parm.Mobility;});
            }
            ExplicitTimeStep(cur_time);
            h_parm->FiniteVolume_Parm.ID_Lambda = 1; //0: General; 1: Special (like all gases or single-phase flows)
            // device_modify([=] AMREX_GPU_DEVICE (){d_parm->FiniteVolume_Parm.ID_Lambda = 1;});
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
            {
                lparm->FiniteVolume_Parm.ID_Lambda = 1;
            });
        }

        if (h_parm->PhaseField_Parm.phase_field){
            if (step >= h_parm->PhaseField_Parm.Switch_Step && t_new[0] >= h_parm->PhaseField_Parm.Switch_Time){
                PhaseField_Evolve();
            }
        }

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
            OutputTime ||
            OutputInterval
            ) {
            last_plot_file_step = step+1;

#ifdef USER_OUTPUT_FUNC 
            if (user_output_int <= 0) {
                UserOutputFunction(*this);
            }
#endif
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

// advance solution to final time
void
Compressible_PhaseField::PhaseField_Evolve ()
{
    Real max_V;
    MaxVelocity(max_V,
                dof_new);

    Real& PF_Lambda = h_parm->PhaseField_Parm.Lambda;
    if (PF_Lambda <= 0.0){
        PF_Lambda = 1.0;
    }

    Real& PF_Eta = h_parm->PhaseField_Parm.Eta;
    if (PF_Eta <= 0.0){
        // The cell size at max_level, which the initial condition also uses, and not
        // that of the finest level at this step, which can be coarser
        const auto dx_Finest = dxFinest();
        PF_Eta = h_parm->PhaseField_Parm.Eta_Multiplier*
                 std::max({ AMREX_D_DECL( dx_Finest[0],
                                          dx_Finest[1],
                                          dx_Finest[2] ) });
    }

    //Real PF_Mobility = 0.5*PF_Eta*max_V;
    Real& PF_Mobility = h_parm->PhaseField_Parm.Mobility;
    if (PF_Mobility == 1.0){
        if (h_parm->PhaseField_Parm.Mechanism == PhaseFieldMechanism::CH){
            PF_Mobility *= 3.0*PF_Eta*PF_Eta*PF_Eta*max_V;
        }
        else{
            PF_Mobility *= 0.5*PF_Eta*max_V;
        }
    }
    else{
        if (h_parm->PhaseField_Parm.Mechanism == PhaseFieldMechanism::CH){
            PF_Mobility *= 0.1*PF_Eta*PF_Eta*PF_Eta;
        }
        else{
            PF_Mobility *= 0.1*PF_Eta;
        }
    }

    const int Num_Grow = NGROW;
    const int Num_Phase = NPHASE;

    //PF boundary conditions
    Vector<BCRec> PF_bcs(Num_Phase);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            //Periodic
            if (geom[0].isPeriodic(idim)){
                PF_bcs[i_Phase].setLo(idim, BCType::int_dir);
                PF_bcs[i_Phase].setHi(idim, BCType::int_dir);
            }
            //Homogeneous Neumann
            else{
                PF_bcs[i_Phase].setLo(idim, BCType::foextrap);
                PF_bcs[i_Phase].setHi(idim, BCType::foextrap);
            }
        }
    }

    //Declare variables
    amrex::Vector<amrex::MultiFab> PFBorder(finest_level+1); //Store Volume Fraction from dof_new
    amrex::Vector<amrex::MultiFab> Div_Jhat(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> Jhat_Flux(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> Grad_QPF(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        PFBorder[lev].define(grids[lev], dmap[lev], Num_Phase, Num_Grow);
        FillPatch(PFBorder[lev], lev, t_new[0],
                  dof_new, t_new,
                  dof_new, t_old,
                  INDEX_VolumeFraction1, 0, Num_Phase,
                  PF_bcs, 0,
                  &lincc_interp);

        Div_Jhat[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        Div_Jhat[lev].setVal(0.0);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            if (h_parm->PhaseField_Parm.ID_ExplicitRC == 1){
                Jhat_Flux[lev][idim].define(ba, dmap[lev], Num_Phase, 0);
                Jhat_Flux[lev][idim].setVal(0.0);
            }
            Grad_QPF[lev][idim].define(ba, dmap[lev], Num_Phase, 0);
            Grad_QPF[lev][idim].setVal(0.0);
        }
    }

    // Compute Div_Jhat and Jhat_Flux
    PhaseFieldMechanism_Evolve(Div_Jhat, Jhat_Flux,
                               PFBorder, PF_bcs,
                               h_parm->PhaseField_Parm);

    // Compute Grad_QPF
    Real PF_Eta_Factor = h_parm->PhaseField_Parm.Eta_Factor;
    if (h_parm->PhaseField_Parm.ID_ExplicitRC == 1){
        if (h_parm->PhaseField_Parm.ID_DegenerateMobility == 1){
            PhaseField_GradQPF(Grad_QPF,
                               Jhat_Flux);
        }
        else{
            PhaseField_GradQPF(Grad_QPF,
                               Jhat_Flux, PFBorder,
                               PF_Eta, PF_Eta_Factor);
        }
    }
    else{
        PhaseField_GradQPF(Grad_QPF,
                           Div_Jhat, PFBorder,
                           PF_Eta, PF_Eta_Factor,
                           h_parm->PhaseField_Parm.ID_Weight,
                           h_parm->LinearSystem_Parm);
    }

    //Reset t_new and told
    for (int lev = 0; lev <= finest_level; lev++){
        t_new[lev] -= dt[0];
        t_old[lev] = t_new[lev];
    }

    //Advance Phase-Field
    Real PF_cur_time = 0.0;
    Real cur_time = t_new[0] + PF_cur_time;

    PF_dt.resize(finest_level+1);

    /* Loop over time steps, stopping when the final time is reached */
    for (int PF_step = 0; PF_step < h_parm->PhaseField_Parm.max_step && PF_cur_time < dt[0]; ++PF_step)
    {
        // amrex::Print() << "\n[PF] Phase Field Coarse STEP " << substep+1 << " starts ..." << std::endl;

        /* Compute the time step on each level */
        PhaseField_ComputeDt(Grad_QPF,PF_cur_time);

        PhaseField_ExplicitTimeStep (cur_time,
                                     Grad_QPF);

        PF_cur_time += PF_dt[0];
        cur_time += PF_dt[0];


//         //=================================================================================
//         //                           Post timestep things
//         //=================================================================================

//         // sum a variable to check conservation and for NaNs 
        Real sum_state = dof_new[0].sum(sum_var);

        // amrex::Print() << "[PF] Coarse SUBSTEP " << substep+1 << " ends." << " PF_Time = " << PF_cur_time
        //                << " PF_DT = " << PF_dt[0] << " Sum(" << vnamesC[sum_var] << ") = " << sum_state << std::endl;

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }


        if (PF_cur_time >= dt[0] - 1.0e-6*dt[0]) {
            amrex::Print() << "[PF] Finished coarse SUBSTEP " << PF_step+1 << "." << " Final PF_Time = " << PF_cur_time
                       << " Final PF_DT = " << PF_dt[0] << " Sum(" << vnamesC[sum_var] << ") = " << sum_state << std::endl;
            break;
        }
        
    }

}


// advance Phase-Field mechanism to one hyperbolic time step
void
Compressible_PhaseField::PhaseFieldMechanism_Evolve (amrex::Vector<amrex::MultiFab>& Div_Jhat,
                                                     amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& Jhat_Flux,
                                                     amrex::Vector<amrex::MultiFab> const& PF,
                                                     amrex::Vector<amrex::BCRec> const& PF_bcs,
                                                     PhaseField_Parameter const& PhaseField_Parm)
{
    const amrex::Real PF_Lambda   = PhaseField_Parm.Lambda;
    const amrex::Real PF_Eta      = PhaseField_Parm.Eta;
    const amrex::Real PF_Mobility = PhaseField_Parm.Mobility;

    //Declare variables
    const int Num_Grow = NGROW;
    const int Num_Phase = NPHASE;
    amrex::Vector<amrex::MultiFab> PF_new(finest_level+1);
    amrex::Vector<amrex::MultiFab> PF_old(finest_level+1);
    amrex::Vector<amrex::MultiFab> PFBorder(finest_level+1);
    amrex::Vector<amrex::MultiFab> Div_J(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> J_Flux(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        PF_new[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        MultiFab::Copy(PF_new[lev], PF[lev], 0, 0, Num_Phase, 0);

        PF_old[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        MultiFab::Copy(PF_old[lev], PF[lev], 0, 0, Num_Phase, 0);

        PFBorder[lev].define(grids[lev], dmap[lev], Num_Phase, Num_Grow);

        Div_J[lev].define(grids[lev], dmap[lev], Num_Phase, 0);

        if (PhaseField_Parm.ID_FluxBasedMechanism == 1){
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
                BoxArray ba = grids[lev];
                ba.surroundingNodes(idim);
                J_Flux[lev][idim].define(ba, dmap[lev], Num_Phase, 0);
            }
        }

        if (PhaseField_Parm.ID_ExplicitRC == 1){
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
                Jhat_Flux[lev][idim].setVal(0.0);
            }
        }
    }

    //Compute PFMechanism_dt
    Real PFMechanism_cfl = 0.9;
    if (PhaseField_Parm.Mechanism == PhaseFieldMechanism::CH){
        PFMechanism_cfl = 0.0000009;
    }
    AMREX_D_TERM(Real dx_min = geom[finest_level].CellSize(0);,
                 Real dy_min = geom[finest_level].CellSize(1);,
                 Real dz_min = geom[finest_level].CellSize(2));
    Real PFMechanism_dt = PFMechanism_cfl/(2.0*PF_Mobility*PF_Lambda)/(AMREX_D_TERM(1.0/dx_min/dx_min,
                                                                                  + 1.0/dy_min/dy_min,
                                                                                  + 1.0/dz_min/dz_min) + 1.0/PF_Eta/PF_Eta);

    //Evolve Phase-Field mechanism
    Real PFMechanism_t = 0.0;
    Real PFMechanism_T = dt[0];
    while (PFMechanism_t < PFMechanism_T){
        //time step
        Real PFMechanism_cur_time = t_new[0] - dt[0] + PFMechanism_t;
        if (PFMechanism_t + PFMechanism_dt > PFMechanism_T){
            PFMechanism_dt = PFMechanism_T - PFMechanism_t;
        }

        //fill PF
        for (int lev = 0; lev <= finest_level; lev++){
            MultiFab::Copy(PF_old[lev], PF_new[lev], 0, 0, Num_Phase, 0);
            FillPatch(PFBorder[lev], lev, PFMechanism_cur_time,
                      PF_new, t_new,
                      PF_new, t_old,
                      0, 0, Num_Phase,
                      PF_bcs, 0,
                      &lincc_interp);
        }
        //Compute Div_J
        if (PhaseField_Parm.ID_FluxBasedMechanism == 1){
            PhaseField_Mechanism(Div_J, J_Flux,
                                 PFBorder, PFBorder,
                                 PhaseField_Parm);
            if (PhaseField_Parm.ID_ExplicitRC == 1){
                for (int lev = 0; lev <= finest_level; ++lev){
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
                        MultiFab::Saxpy(Jhat_Flux[lev][idim], PFMechanism_dt, J_Flux[lev][idim], 0, 0, Num_Phase, 0);
                    }
                }
            }
        }
        else{
            PhaseField_Mechanism(Div_J,
                                 PFBorder, PFBorder,
                                 PhaseField_Parm);
        }
        //Update PF
        for (int lev = 0; lev <= finest_level; lev++){
            MultiFab::LinComb(PF_new[lev], Real(1.0), PF_old[lev], 0, PFMechanism_dt, Div_J[lev], 0, 0, Num_Phase, 0);
        }

        PFMechanism_t += PFMechanism_dt;
    }
    for (int lev = 0; lev <= finest_level; lev++){
        if (PhaseField_Parm.ID_ExplicitRC == 1){
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
                Jhat_Flux[lev][idim].mult((1.0/dt[0]), 0, Num_Phase, 0);
            }
        }
        else{
            MultiFab::LinComb(Div_Jhat[lev], (1.0/dt[0]), PF_new[lev], 0, (-1.0/dt[0]), PF[lev], 0, 0, Num_Phase, 0);
        }
    }

}


//Advance state vector at a level with different RK schemes
void 
Compressible_PhaseField::PhaseField_AdvanceAtLevel (Vector<MultiFab>& mf_new, 
                                                    Vector<MultiFab>& mf_old, 
                                                    Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                                    int lev, 
                                                    Real time, 
                                                    Real dt_lev, 
                                                    int ncycle)
{

    if (ID_TimeIntegrator == "ForwardEuler"){
        PhaseField_ForwardEuler(mf_new,
                            mf_old, 
                            Grad_Q, 
                            lev,
                            time,
                            dt_lev,
                            ncycle);

    }  else if (ID_TimeIntegrator == "TVD-RK2"){
        PhaseField_SecondOrderSSPRK(mf_new,
                                    mf_old,
                                    Grad_Q, 
                                    lev,
                                    time,
                                    dt_lev,
                                    ncycle);

    }  else if (ID_TimeIntegrator == "TVD-RK3"){
        PhaseField_ThirdOrderSSPRK(mf_new,
                                    mf_old,
                                    Grad_Q, 
                                    lev,
                                    time,
                                    dt_lev,
                                    ncycle);

    } else if (ID_TimeIntegrator == "RK4"){
        PhaseField_FourthOrderRK(mf_new,
                                    mf_old,
                                    Grad_Q, 
                                    lev,
                                    time,
                                    dt_lev,
                                    ncycle);
    }

}


// PhaseField_AdvanceAllLevels (the step without subcycling): Class_Compressible_PhaseField_StageCoupling.cpp


// Phase-Field face fluxes of level lev from the state U with ghost cells and the face gradients
// Grad_Q of the level: in fluxes the conservative fluxes per unit area, in fluxes_nc the face
// quantities of the Phase-Field non-conservative terms, and the speeds in PF_c_max[lev] (reset at
// stage 0). As in ComputeFaceFluxes each tile computes the faces of nodaltilebox only
void
Compressible_PhaseField::PhaseField_ComputeFaceFluxes (MultiFab const& U,
                                                       Array<MultiFab,AMREX_SPACEDIM> const& Grad_Q,
                                                       Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                                       Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                                       int lev,
                                                       int stage)
{
BL_PROFILE("PhaseField_ComputeFaceFluxes()");

    MultiFab& c_max_lev = PF_c_max[lev];
    if (stage == 0){
        c_max_lev.setVal(0.0);
    }

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
    amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
    amrex::ignore_unused(dy, dz, fluxes_nc);

    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(U,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Pull the data into an array
            Array4<Real const> statein   = U.const_array(mfi);
            Array4<Real      > c_max_new = c_max_lev.array(mfi);
{BL_PROFILE("PhaseField_compute_dUdt_FV::{ computing the Phase-Field fluxes }");
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                Array4<Real const> GradX_Q = Grad_Q[idim].const_array(mfi);
                Array4<Real      > flux    = fluxes[idim].array(mfi);
#if (NONCONSERVATIVE == true)
                Array4<Real      > fluxNC  = fluxes_nc[idim].array(mfi);
#endif

                // the faces of this tile only (the high face of a box with its last tile)
                const Box bx_face = mfi.nodaltilebox(idim);
                amrex::ParallelFor(bx_face,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Array<Real,NPHASE> grad_Q;
                    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
                        grad_Q[i_Phase] = GradX_Q(i,j,k,i_Phase);
                    }

                    Array<Real,NSTATE> fhat;
                    Array<Real,NC_TERMS> fhatNC;
                    Real cmax;

                    FVM_Conservative2FluxPhaseField_K(i, j, k,
                                                      fhat, fhatNC, cmax,
                                                      statein, grad_Q, *lparm,
                                                      dX[idim], idim,
                                                      lparm->FiniteVolume_Parm, lparm->PhaseField_Parm);

                    c_max_new(i,j,k,0) = std::max(cmax, c_max_new(i,j,k,0));
                    for (int iState = 0; iState < NSTATE; iState++){
                        flux(i,j,k,iState) = -fhat[iState];
                    }
#if (NONCONSERVATIVE == true)
                    for (int iState = 0; iState < NC_TERMS; iState++){
                        fluxNC(i,j,k,iState) = fhatNC[iState];
                    }
#endif
                });
            } // end idim
}
        } // end mfi
    } // end omp
}


// Phase-Field right-hand side dUdt of level lev (overwritten) from the face fluxes of
// PhaseField_ComputeFaceFluxes: the conservative surface integral and the Phase-Field
// non-conservative terms, from the same face MultiFabs
void
Compressible_PhaseField::PhaseField_FluxDivergence (MultiFab& dUdt_mf,
                                                    MultiFab const& U,
                                                    Array<MultiFab,AMREX_SPACEDIM> const& fluxes,
                                                    Array<MultiFab,AMREX_SPACEDIM> const& fluxes_nc,
                                                    int lev)
{
BL_PROFILE("PhaseField_FluxDivergence()");

    auto const prob_lo = Geom(lev).ProbLoArray();
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
    amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
    amrex::ignore_unused(dy, dz, fluxes_nc);

    Parm const* lparm = d_parm;

    const int lID_COORDSYS = ID_COORDSYS;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(dUdt_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            Array4<Real const> statein = U.const_array(mfi);
            Array4<Real      > dUdt    = dUdt_mf.array(mfi);

            AMREX_D_TERM(Array4<Real const> fluxx_c = fluxes[0].const_array(mfi);,
                         Array4<Real const> fluxy_c = fluxes[1].const_array(mfi);,
                         Array4<Real const> fluxz_c = fluxes[2].const_array(mfi));
{BL_PROFILE("COMPAS::PhaseField_compute_dUdt_FV::surface_integral()");
            // Do a conservative update
            // Forward Euler
            // ===========================
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (lID_COORDSYS == 0){
                    for (int iState = 0; iState < NSTATE; iState++){
                        FVM_SurfaceIntegral(i, j, k, iState,
                                         dUdt, statein,
                            AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
                            AMREX_D_DECL(dx,dy,dz));
                    }
                }
                else {
                    FVM_SurfaceIntegral_Coord(i,j,k,
                                           dUdt,statein,
                              AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
                                           lID_COORDSYS,
                                           dX,
                                           prob_lo);
                }
            });
}

#if (NONCONSERVATIVE == true)
            AMREX_D_TERM(Array4<Real const> fluxxNC_c = fluxes_nc[0].const_array(mfi);,
                         Array4<Real const> fluxyNC_c = fluxes_nc[1].const_array(mfi);,
                         Array4<Real const> fluxzNC_c = fluxes_nc[2].const_array(mfi));
{BL_PROFILE("COMPAS::PhaseField_compute_dUdt_FV::PhaseField_surface_integral_nc()");
            // Add the non-conservative contribution
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                PhaseField_FVM_SurfaceIntegral_NC(i, j, k,
                                               dUdt, statein,
                                  AMREX_D_DECL(fluxxNC_c, fluxyNC_c, fluxzNC_c),
                                  AMREX_D_DECL(dx       , dy       , dz       ),
                                               *lparm);
            });
}
#endif
        } // end mfi
    } // end omp
}


//Advance state vector at a level with the Phase-Field flux vector (per-level integrators of the
//subcycled step; with do_reflux the area-weighted fluxes go to the flux registers)
void
Compressible_PhaseField::PhaseField_compute_dUdt_FV (MultiFab& mf_new,
                                                     MultiFab& mf_old,
                                                     Array<MultiFab,AMREX_SPACEDIM> const& Grad_Q,
                                                     int lev,
                                                     Real dt_lev,
                                                     int ncycle,
                                                     int stage,
                                                     FluxRegister* fr_as_crse,
                                                     FluxRegister* fr_as_fine,
                                                     FluxRegister* fr_as_crse_nc,
                                                     FluxRegister* fr_as_fine_nc)
{
BL_PROFILE("PhaseField_compute_dUdt_FV()");
    amrex::ignore_unused(ncycle, fr_as_crse_nc, fr_as_fine_nc);

    // construct the fluxes
    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    Array<MultiFab,AMREX_SPACEDIM> fluxes_nc;
    DefineFaceFluxes(fluxes, fluxes_nc, lev);

    PhaseField_ComputeFaceFluxes(mf_old, Grad_Q, fluxes, fluxes_nc, lev, stage);

    PhaseField_FluxDivergence(mf_new, mf_old, fluxes, fluxes_nc, lev);

    if (do_reflux)
    {
        ScaleFaceFluxesByArea(fluxes, fluxes_nc, lev);
    }

    if (fr_as_crse) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = -dt_lev;
            fr_as_crse->CrseInit(fluxes[idim], idim, 0, 0, NSTATE, scale, FluxRegister::ADD);
#if (NONCONSERVATIVE == true)
            // fr_as_crse_nc->CrseInit(fluxes_nc[idim], idim, 0, 0, NC_FLUX_COMP, scale, FluxRegister::ADD);
            fr_as_crse_nc->CrseInit(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale, FluxRegister::ADD);
#endif
        }
    }

    if (fr_as_fine) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = dt_lev;
            fr_as_fine->FineAdd(fluxes[idim], idim, 0, 0, NSTATE, scale);
#if (NONCONSERVATIVE == true)
            // fr_as_fine_nc->FineAdd(fluxes_nc[idim], idim, 0, 0, NC_FLUX_COMP, scale);
            fr_as_fine_nc->FineAdd(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale);
#endif
        }
    }
}

// Compute Div_J from the conservative Allen-Cahn model
void
Compressible_PhaseField::PhaseField_Mechanism(amrex::Vector<amrex::MultiFab>& Div_J,
                                              amrex::Vector<amrex::MultiFab> const& PFBorder,
                                              amrex::Vector<amrex::MultiFab> const& PsiBorder,
                                              PhaseField_Parameter const& PhaseField_Parm)
{
    const int Num_Phase    = PFBorder[0].nComp();
    const int ID_Mechanism = PhaseField_Parm.ID_Mechanism;
    const int ID_Potential = PhaseField_Parm.ID_Potential;
    const int ID_Weight    = PhaseField_Parm.ID_Weight;

    const amrex::Real Mobility = PhaseField_Parm.Mobility;
    const amrex::Real Lambda   = PhaseField_Parm.Lambda;
    const amrex::Real Eta      = PhaseField_Parm.Eta;

    //Define Mechanism_PF
    Vector<MultiFab> Mechanism_PF(finest_level+1);
    //Define Weight_PF
    Vector<MultiFab> Weight_PF(finest_level+1);
    //Define Mobility_PF
    Vector<MultiFab> Mobility_PF(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        Mechanism_PF[lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        Weight_PF   [lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        Mobility_PF [lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        Mobility_PF [lev].setVal(Mobility);
    }
    //Define NormalBorder
    Vector<MultiFab> NormalBorder(finest_level+1);

    // =======================================================
    // Compute Div_J for all phases at all levels
    // =======================================================
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        // =======================================================
        //Compute Normal vector at cell centers at all levels from PF
        // =======================================================
        PhaseField_Normal(NormalBorder,
                          PsiBorder,
                          i_Phase);

        // =======================================================
        // Compute chemical potential and weight at all levels
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(Div_J[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx  = mfi.tilebox();

                    Array4<Real const> PFArray        = PFBorder    [lev].const_array(mfi);
                    Array4<Real const> NormalArray    = NormalBorder[lev].const_array(mfi);
                    Array4<Real const> MobilityArray  = Mobility_PF [lev].const_array(mfi);
                    Array4<Real      > MechanismArray = Mechanism_PF[lev].array(mfi);
                    Array4<Real      > WeightArray    = Weight_PF   [lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        amrex::Real PF_Mechanism;
                        if (ID_Mechanism == 0){
                            FDM_PhaseFieldMechanism_ACVanDerWaals_K(i, j, k, i_Phase,
                                                                    PF_Mechanism,
                                                                    PFArray,
                                                                    MobilityArray, Lambda, Eta,
                                                       AMREX_D_DECL(dx, dy, dz),
                                                                    ID_Potential);
                        }
                        else if (ID_Mechanism == 1){
                            FDM_PhaseFieldMechanism_ACAdv_K(i, j, k, i_Phase,
                                                            PF_Mechanism,
                                                            PFArray, NormalArray,
                                                            MobilityArray, Lambda, Eta,
                                               AMREX_D_DECL(dx, dy, dz),
                                                            ID_Potential);
                        }
                        else if (ID_Mechanism == 2){
                            FDM_PhaseFieldMechanism_ACAdv_NormalFace_K(i, j, k, i_Phase,
                                                                       PF_Mechanism,
                                                                       PFArray,
                                                                       MobilityArray, Lambda, Eta,
                                                          AMREX_D_DECL(dx, dy, dz),
                                                                       ID_Potential);
                        }
                        MechanismArray(i,j,k,0) = PF_Mechanism;

                        //Weight
                        int ID_Derivative = 0;
                        amrex::Real PF_Weight;
                        WeightFunction_PhaseField(PF_Weight,
                                                  PFArray(i,j,k,i_Phase),
                                                  ID_Derivative, ID_Weight);
                        WeightArray(i,j,k,0) = PF_Weight;
                    });
                } // end mfi
            } // end omp
        } // end lev

        // ==========================================================
        // Compute domain integral of ChemicalPotential and Weight_PF
        // ==========================================================
        AverageDown(Mechanism_PF);
        AverageDown(Weight_PF);
        amrex::Real Integral_Mechanism = Mechanism_PF[0].sum(0);
        amrex::Real Integral_Weight    = Weight_PF[0].sum(0);

        // ==========================================================
        // Compute Div_J at all levels
        // ==========================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(Div_J[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx  = mfi.tilebox();

                    Array4<Real const> MechanismArray = Mechanism_PF[lev].const_array(mfi);
                    Array4<Real const> WeightArray    = Weight_PF   [lev].const_array(mfi);
                    Array4<Real      > Div_JArray     = Div_J       [lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        amrex::Real PF_Mechanism_AC = MechanismArray(i,j,k,0);
                        amrex::Real PF_Weight       = WeightArray(i,j,k,0);
                        amrex::Real PF_Mechanism;
                        PhaseFieldMechanism_CAC(PF_Mechanism,
                                                PF_Mechanism_AC, Integral_Mechanism,
                                                PF_Weight      , Integral_Weight);
                        Div_JArray(i,j,k,i_Phase) = PF_Mechanism;
                    });
                } // end mfi
            } // end omp
        } // end lev
    } // end i_Phase
}


void
Compressible_PhaseField::PhaseField_Mechanism(amrex::Vector<amrex::MultiFab>& Div_J,
                                              amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& J_Flux,
                                              amrex::Vector<amrex::MultiFab> const& PFBorder,
                                              amrex::Vector<amrex::MultiFab> const& PsiBorder,
                                              PhaseField_Parameter const& PhaseField_Parm)
{
    const int Num_Phase    = PFBorder[0].nComp();
    const int ID_Mechanism = PhaseField_Parm.ID_Mechanism;
    const int ID_Potential = PhaseField_Parm.ID_Potential;

    const amrex::Real Mobility = PhaseField_Parm.Mobility;
    const amrex::Real Lambda   = PhaseField_Parm.Lambda;
    const amrex::Real Eta      = PhaseField_Parm.Eta;

    const PhaseFieldMechanism Mechanism = PhaseField_Parm.Mechanism;
    if (Mechanism == PhaseFieldMechanism::CDI){
        PhaseField_CDI(J_Flux,
                       PFBorder, PsiBorder,
                       Mobility, Eta,
                       ID_Mechanism);
    }
    else if (Mechanism == PhaseFieldMechanism::ACDI){
        PhaseField_ACDI(J_Flux,
                        PFBorder, PsiBorder,
                        Mobility, Eta,
                        ID_Mechanism);
    }
    else if (Mechanism == PhaseFieldMechanism::CH){
        PhaseField_CH(J_Flux,
                      PFBorder, PsiBorder,
                      Mobility, Lambda, Eta,
                      ID_Mechanism, ID_Potential);
    }
    AverageDownFaces(J_Flux);

    // =======================================================
    // Compute Div_J for all phases at all levels
    // =======================================================
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        // ==========================================================
        // Compute Div_J at all levels
        // ==========================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(Div_J[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();

       AMREX_D_TERM(Array4<Real const> J_FluxArray_x = J_Flux[lev][0].const_array(mfi);,
                    Array4<Real const> J_FluxArray_y = J_Flux[lev][1].const_array(mfi);,
                    Array4<Real const> J_FluxArray_z = J_Flux[lev][2].const_array(mfi));
                    Array4<Real      > Div_JArray    = Div_J [lev].array(mfi);

                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        amrex::Real Div_f;
                        FVM_Divergence_K(i, j, k, i_Phase,
                                         Div_f,
                            AMREX_D_DECL(J_FluxArray_x, J_FluxArray_y, J_FluxArray_z),
                            AMREX_D_DECL(           dx,            dy,            dz));
                        Div_JArray(i,j,k,i_Phase) = Div_f;
                    });
                } // end mfi
            } // end omp
        } // end lev
    } // end i_Phase
}


// Compute J_Flux from the Cahn-Hilliard model
void
Compressible_PhaseField::PhaseField_CH(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& J_Flux,
                                       amrex::Vector<amrex::MultiFab> const& PFBorder,
                                       amrex::Vector<amrex::MultiFab> const& PsiBorder,
                                       amrex::Real const& Mobility,
                                       amrex::Real const& Lambda,
                                       amrex::Real const& Eta,
                                       int ID_Mechanism,
                                       int ID_Potential)
{
    const int Num_Phase = PFBorder[0].nComp();

    //Define ChemicalPotential
    Vector<MultiFab> ChemicalPotential(finest_level+1);
    //Define Mobility
    Vector<MultiFab> Mobility_PF(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        ChemicalPotential[lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        Mobility_PF      [lev] = MultiFab(grids[lev], dmap[lev], 1, 1);
        Mobility_PF      [lev].setVal(Mobility);
    }
    //Define ChemicalPotential with ghost cells
    MultiFab ChemicalPotentialBorder;
    //ChemicalPotential boundary conditions
    Vector<BCRec> ChemicalPotential_bcs(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        //Periodic
        if (geom[0].isPeriodic(idim)){
            ChemicalPotential_bcs[0].setLo(idim, BCType::int_dir);
            ChemicalPotential_bcs[0].setHi(idim, BCType::int_dir);
        }
        //Homogeneous Neumann
        else{
            ChemicalPotential_bcs[0].setLo(idim, BCType::foextrap);
            ChemicalPotential_bcs[0].setHi(idim, BCType::foextrap);
        }
    }
    //Define NormalBorder
    Vector<MultiFab> NormalBorder(finest_level+1);

    // =======================================================
    // Compute J_Flux for all phases at all levels
    // =======================================================
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        // =======================================================
        //Compute Normal vector at cell centers at all levels from PF
        // =======================================================
        PhaseField_Normal(NormalBorder,
                          PsiBorder,
                          i_Phase);

        // =======================================================
        // Compute ChemicalPotential at all levels
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(ChemicalPotential[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx  = mfi.tilebox();
                
                    Array4<Real const> PFArray                = PFBorder         [lev].const_array(mfi);
                    Array4<Real const> NormalArray            = NormalBorder     [lev].const_array(mfi);
                    Array4<Real      > ChemicalPotentialArray = ChemicalPotential[lev].array(mfi);
                    
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        amrex::Real ChemicalPotential_c;
                        if (ID_Mechanism == 5){
                            FDM_ChemicalPotential_VanDerWaals_K(i, j, k, i_Phase,
                                                                ChemicalPotential_c,
                                                                PFArray,
                                                                Lambda, Eta,
                                                   AMREX_D_DECL(dx, dy, dz),
                                                                ID_Potential);
                        }
                        else if (ID_Mechanism == 6){
                            FDM_ChemicalPotential_Adv_K(i, j, k, i_Phase,
                                                        ChemicalPotential_c,
                                                        PFArray, NormalArray,
                                                        Lambda, Eta,
                                           AMREX_D_DECL(dx, dy, dz),
                                                        ID_Potential);
                        }
                        else if (ID_Mechanism == 7){
                            FDM_ChemicalPotential_Adv_NormalFace_K(i, j, k, i_Phase,
                                                                   ChemicalPotential_c,
                                                                   PFArray,
                                                                   Lambda, Eta,
                                                      AMREX_D_DECL(dx, dy, dz),
                                                                   ID_Potential);
                        }
                        ChemicalPotentialArray(i,j,k,0) = ChemicalPotential_c;
                    });
                } // end mfi
            } // end omp
        } // end lev
        
        // =======================================================
        // Compute J_Flux at all levels
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            ChemicalPotentialBorder = MultiFab(grids[lev], dmap[lev], 1, 1);
            FillPatch(ChemicalPotentialBorder, lev, t_new[0],
                      ChemicalPotential, t_new,
                      ChemicalPotential, t_new,
                      0, 0, 1,
                      ChemicalPotential_bcs, 0,
                      &lincc_interp);

            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(ChemicalPotentialBorder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();

                    Array4<Real const> ChemicalPotentialArray = ChemicalPotentialBorder.const_array(mfi);
                    Array4<Real const> MobilityArray          = Mobility_PF       [lev].const_array(mfi);

                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        Array4<Real> J_FluxArray = J_Flux[lev][idim].array(mfi);
                    
                        // the faces of this tile only, so that no two OpenMP threads write one face
                        amrex::ignore_unused(bx);
                        const Box bx_face = mfi.nodaltilebox(idim);
                        amrex::ParallelFor(bx_face,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            amrex::Real Mobility_Face;
                            Interpolation2Face_Linear_K(i, j, k, 0,
                                                        Mobility_Face,
                                                        MobilityArray,
                                                        idim);

                            amrex::Real Grad_ChemicalPotential;
                            FDM_Gradient_FaceFromCenter_K(i, j, k, 0,
                                                          Grad_ChemicalPotential,
                                                          ChemicalPotentialArray,
                                                          dX[idim],
                                                          idim);
                            J_FluxArray(i,j,k,i_Phase) = Mobility_Face*Grad_ChemicalPotential;
                        });
                    } // end idim
                } // end mfi
            } // end omp
        } // end lev
    } // end i_Phase
}


// Compute J_Flux from the conservative Diffuse-Interface (Phase-Field) model
void
Compressible_PhaseField::PhaseField_CDI(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& J_Flux,
                                        amrex::Vector<amrex::MultiFab> const& PFBorder,
                                        amrex::Vector<amrex::MultiFab> const& PsiBorder,
                                        amrex::Real const& Mobility,
                                        amrex::Real const& Eta,
                                        int ID_Mechanism)
{
    const int Num_Phase = PFBorder[0].nComp();

    amrex::Vector<amrex::MultiFab> Mobility_PF(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        Mobility_PF[lev].define(grids[lev], dmap[lev], 1, 1);
        Mobility_PF[lev].setVal(Mobility);
    }

    amrex::Vector<amrex::MultiFab> NormalBorder(finest_level+1);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        // =======================================================
        // Compute NormalBorder at all levels with PsiBorder
        // =======================================================
        PhaseField_Normal(NormalBorder,
                          PsiBorder,
                          i_Phase);

        // =======================================================
        // Compute J_Flux at all levels
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(PFBorder[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();

                    Array4<Real const> PFArray       = PFBorder    [lev].const_array(mfi);
                    Array4<Real const> NormalArray   = NormalBorder[lev].const_array(mfi);
                    Array4<Real const> MobilityArray = Mobility_PF [lev].const_array(mfi);

                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        Array4<Real> J_FluxArray = J_Flux[lev][idim].array(mfi);
                    
                        // the faces of this tile only, so that no two OpenMP threads write one face
                        amrex::ignore_unused(bx);
                        const Box bx_face = mfi.nodaltilebox(idim);
                        amrex::ParallelFor(bx_face,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            amrex::Real PF_Flux;
                            if (ID_Mechanism == 0){
                                FDM_PhaseFieldFlux_CDI_K(i, j, k, i_Phase,
                                                         PF_Flux,
                                                         PFArray, NormalArray,
                                                         MobilityArray, Eta,
                                                         dX[idim],
                                                         idim);
                            }
                            else if (ID_Mechanism == 1){
                                FDM_PhaseFieldFlux_CDI_CompressionInterpolation_K(i, j, k, i_Phase,
                                                                                  PF_Flux,
                                                                                  PFArray, NormalArray,
                                                                                  MobilityArray, Eta,
                                                                                  dX[idim],
                                                                                  idim);
                            }
                            else if (ID_Mechanism == 2){
                                FDM_PhaseFieldFlux_CDI_NormalFace_K(i, j, k, i_Phase,
                                                                    PF_Flux,
                                                                    PFArray,
                                                                    MobilityArray, Eta,
                                                       AMREX_D_DECL(dx, dy, dz),
                                                                    idim);
                            }
                            J_FluxArray(i,j,k,i_Phase) = PF_Flux;
                        });
                    } // end idim
                } // end mfi
            } // end omp
        } // end lev
    } // end i_Phase
}


// Compute J_Flux from the accurate conservative Diffuse-Interface model
void
Compressible_PhaseField::PhaseField_ACDI(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& J_Flux,
                                         amrex::Vector<amrex::MultiFab> const& PFBorder,
                                         amrex::Vector<amrex::MultiFab> const& PsiBorder,
                                         amrex::Real const& Mobility,
                                         amrex::Real const& Eta,
                                         int ID_Mechanism)
{
    const int Num_Phase = PFBorder[0].nComp();

    //Define Mobility
    amrex::Vector<amrex::MultiFab> Mobility_PF(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        Mobility_PF[lev].define(grids[lev], dmap[lev], 1, 1);
        Mobility_PF[lev].setVal(Mobility);
    }

    //Define LS
    amrex::Vector<amrex::MultiFab> LS(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        LS[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
    }
    //LS boundary conditions
    Vector<BCRec> LS_bcs(Num_Phase);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            //Periodic
            if (geom[0].isPeriodic(idim)){
                LS_bcs[i_Phase].setLo(idim, BCType::int_dir);
                LS_bcs[i_Phase].setHi(idim, BCType::int_dir);
            }
            //Homogeneous Neumann
            else{
                LS_bcs[i_Phase].setLo(idim, BCType::foextrap);
                LS_bcs[i_Phase].setHi(idim, BCType::foextrap);
            }
        }
    }

    // =======================================================
    // Compute LS at all level
    // =======================================================
    for (int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(LS[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx  = mfi.tilebox();
                
                Array4<Real const> PFArray  = PFBorder[lev].const_array(mfi);
                Array4<Real      > LSArray  = LS      [lev].array(mfi);
                
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
                        amrex::Real PF_c = PFArray(i,j,k,i_Phase);

                        int ID_Derivative = 0;
                        amrex::Real LS_c;
                        Mapping_Backward(LS_c,
                                         PF_c, Eta,
                                         ID_Derivative);
                        LSArray(i,j,k,i_Phase) = LS_c;
                    } // end i_Phase
                });
            } // end mfi
        } // end omp
    } // end lev
    // =======================================================
    // Fill LSBorder
    // =======================================================
    amrex::Vector<amrex::MultiFab> LSBorder(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        LSBorder[lev].define(grids[lev], dmap[lev], Num_Phase, 1);
        FillPatch(LSBorder[lev], lev, t_new[0],
                  LS, t_new,
                  LS, t_new,
                  0, 0, Num_Phase,
                  LS_bcs, 0,
                  &lincc_interp);
    }

    amrex::Vector<amrex::MultiFab> NormalBorder(finest_level+1);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        // =======================================================
        // Compute NormalBorder at all levels with PsiBorder
        // =======================================================
        PhaseField_Normal(NormalBorder,
                          LSBorder,
                          i_Phase);

        // =======================================================
        // Compute J_Flux at all levels
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(PFBorder[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();

                    Array4<Real const> PFArray       = PFBorder    [lev].const_array(mfi);
                    Array4<Real const> NormalArray   = NormalBorder[lev].const_array(mfi);
                    Array4<Real const> MobilityArray = Mobility_PF [lev].const_array(mfi);

                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        Array4<Real> J_FluxArray = J_Flux[lev][idim].array(mfi);
                    
                        // the faces of this tile only, so that no two OpenMP threads write one face
                        amrex::ignore_unused(bx);
                        const Box bx_face = mfi.nodaltilebox(idim);
                        amrex::ParallelFor(bx_face,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            amrex::Real PF_Flux;
                            if (ID_Mechanism == 3){
                                FDM_PhaseFieldFlux_ACDI_K(i, j, k, i_Phase,
                                                          PF_Flux,
                                                          PFArray, NormalArray,
                                                          MobilityArray, Eta,
                                                          dX[idim],
                                                          idim);
                            }
                            else if (ID_Mechanism == 4){
                                FDM_PhaseFieldFlux_ACDI_NormalFace_K(i, j, k, i_Phase,
                                                                     PF_Flux,
                                                                     PFArray,
                                                                     MobilityArray, Eta,
                                                        AMREX_D_DECL(dx, dy, dz),
                                                                     idim);
                            }
                            J_FluxArray(i,j,k,i_Phase) = PF_Flux;
                        });
                    } // end idim
                } // end mfi
            } // end omp
        } // end lev
    } // end i_Phase
}


// Compute diffused Heaviside function based on the Phase-Field function
void
Compressible_PhaseField::PhaseField_Diffusion(amrex::Vector<amrex::MultiFab>& HeavisideBorder,
                                              amrex::Vector<amrex::MultiFab> const& PF,
                                              amrex::Real const& Eta)
{
    const int Num_Phase = PF[0].nComp();
    //Define Heaviside
    amrex::Vector<amrex::MultiFab> Heaviside(finest_level+1);
    amrex::Vector<amrex::MultiFab> dHeavisidedt(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        Heaviside[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        dHeavisidedt[lev].define(grids[lev], dmap[lev], Num_Phase, 0);

        HeavisideBorder[lev].define(grids[lev], dmap[lev], Num_Phase, 1);
    }
    //Heaviside boundary conditions
    Vector<BCRec> Heaviside_bcs(Num_Phase);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            //Periodic
            if (geom[0].isPeriodic(idim)){
                Heaviside_bcs[i_Phase].setLo(idim, BCType::int_dir);
                Heaviside_bcs[i_Phase].setHi(idim, BCType::int_dir);
            }
            //Homogeneous Neumann
            else{
                Heaviside_bcs[i_Phase].setLo(idim, BCType::foextrap);
                Heaviside_bcs[i_Phase].setHi(idim, BCType::foextrap);
            }
        }
    }

    // =======================================================
    // Compute Heaviside at all level
    // =======================================================
    amrex::Real Heaviside_max = 1.0;
    amrex::Real Heaviside_min = 0.0;
    amrex::Real Heaviside_mid = 0.5;
    for (int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(Heaviside[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx  = mfi.tilebox();
                
                Array4<Real const> PFArray        = PF       [lev].const_array(mfi);
                Array4<Real      > HeavisideArray = Heaviside[lev].array(mfi);
                
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
                        amrex::Real PF_c = PFArray(i,j,k,i_Phase);
                        amrex::Real Heaviside_c;
                        if (PF_c > Heaviside_mid){
                            Heaviside_c = Heaviside_max;
                        }
                        else if (PF_c < Heaviside_mid){
                            Heaviside_c = Heaviside_min;
                        }
                        else{
                            Heaviside_c = Heaviside_mid;
                        }
                        HeavisideArray(i,j,k,i_Phase) = Heaviside_c;
                    } // end i_Phase
                });
            } // end mfi
        } // end omp
    } // end lev

    // =======================================================
    // Diffuse Heaviside at all level
    // =======================================================
    amrex::Real Diff_L = 7.0*Eta;
    amrex::Real Diff_cfl = 0.99;
    amrex::Real Diff = 1.0;
    amrex::Real Diff_T = Diff_L*Diff_L/Diff;
    AMREX_D_TERM(amrex::Real dx_min = geom[finest_level].CellSize(0);,
                 amrex::Real dy_min = geom[finest_level].CellSize(1);,
                 amrex::Real dz_min = geom[finest_level].CellSize(2));
    amrex::Real Diff_dt = Diff_cfl/(2.0*Diff)/(AMREX_D_TERM(1.0/dx_min/dx_min,
                                                          + 1.0/dy_min/dy_min,
                                                          + 1.0/dz_min/dz_min));
    int Diff_Step = int(Diff_T/Diff_dt + 1.0);
    for (int Diff_step = 0; Diff_step <= Diff_Step; Diff_step++){
        //Fill HeavisideBorder
        for (int lev = 0; lev <= finest_level; lev++){
            FillPatch(HeavisideBorder[lev], lev, t_new[0],
                      Heaviside, t_new,
                      Heaviside, t_new,
                      0, 0, Num_Phase,
                      Heaviside_bcs, 0,
                      &lincc_interp);
        }
        //Compute dHeavisidedt
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(dHeavisidedt[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx  = mfi.tilebox();
                
                    Array4<Real const> HeavisideArray    = HeavisideBorder[lev].const_array(mfi);
                    Array4<Real      > dHeavisidedtArray = dHeavisidedt   [lev].array(mfi);
                
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
                            amrex::Real Laplace_Heaviside;
                            FDM_Laplace_K(i, j, k, i_Phase,
                                          Laplace_Heaviside,
                                          HeavisideArray,
                             AMREX_D_DECL(dx, dy, dz));
                            dHeavisidedtArray(i,j,k,i_Phase) = Diff*Laplace_Heaviside;
                        } // end i_Phase
                    });
                } // end mfi
            } // end omp
        } // end lev
        //Update Heaviside
        for (int lev = 0; lev <= finest_level; lev++){
            amrex::MultiFab::Saxpy(Heaviside[lev], Diff_dt, dHeavisidedt[lev], 0, 0, Num_Phase, 0);
        }
    } // end Diff_step
    //Fill HeavisideBorder
    for (int lev = 0; lev <= finest_level; lev++){
        FillPatch(HeavisideBorder[lev], lev, t_new[0],
                  Heaviside, t_new,
                  Heaviside, t_new,
                  0, 0, Num_Phase,
                  Heaviside_bcs, 0,
                  &lincc_interp);
    }
}


// Compute signed distance function based on the Phase-Field function
void
Compressible_PhaseField::PhaseField_Reinitialization(amrex::Vector<amrex::MultiFab>& LSBorder,
                                                     amrex::Vector<amrex::MultiFab> const& PF,
                                                     amrex::Real const& Eta)
{
    const int Num_Phase = PF[0].nComp();
    //Define LS
    amrex::Vector<amrex::MultiFab> LS0(finest_level+1);
    amrex::Vector<amrex::MultiFab> LS(finest_level+1);
    amrex::Vector<amrex::MultiFab> dLSdt(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        LS0[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        LS[lev].define(grids[lev], dmap[lev], Num_Phase, 0);
        dLSdt[lev].define(grids[lev], dmap[lev], Num_Phase, 0);

        LSBorder[lev].define(grids[lev], dmap[lev], Num_Phase, 1);
    }
    //LS boundary conditions
    Vector<BCRec> LS_bcs(Num_Phase);
    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            //Periodic
            if (geom[0].isPeriodic(idim)){
                LS_bcs[i_Phase].setLo(idim, BCType::int_dir);
                LS_bcs[i_Phase].setHi(idim, BCType::int_dir);
            }
            //Homogeneous Neumann
            else{
                LS_bcs[i_Phase].setLo(idim, BCType::foextrap);
                LS_bcs[i_Phase].setHi(idim, BCType::foextrap);
            }
        }
    }

    // =======================================================
    // Compute LS at all level
    // =======================================================
    for (int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(LS[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx  = mfi.tilebox();
                
                Array4<Real const> PFArray  = PF [lev].const_array(mfi);
                Array4<Real      > LS0Array = LS0[lev].array(mfi);
                Array4<Real      > LSArray  = LS [lev].array(mfi);
                
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
                        amrex::Real PF_c = PFArray(i,j,k,i_Phase);

                        int ID_Derivative = 0;
                        amrex::Real LS_c;
                        Mapping_Backward(LS_c,
                                         PF_c, Eta,
                                         ID_Derivative);
                        LS0Array(i,j,k,i_Phase) = LS_c;
                        LSArray(i,j,k,i_Phase) = LS_c;
                    } // end i_Phase
                });
            } // end mfi
        } // end omp
    } // end lev

    // =======================================================
    // Reinitialization at all level
    // =======================================================
    amrex::Real Rein_cfl = 0.99;
    AMREX_D_TERM(amrex::Real dx_min = geom[finest_level].CellSize(0);,
                 amrex::Real dy_min = geom[finest_level].CellSize(1);,
                 amrex::Real dz_min = geom[finest_level].CellSize(2));
    amrex::Real Rein_Diff = std::max({AMREX_D_DECL(dx_min, dy_min, dz_min)});
    amrex::Real Rein_dt = Rein_cfl/(2.0*Rein_Diff)/(AMREX_D_TERM(1.0/dx_min/dx_min,
                                                               + 1.0/dy_min/dy_min,
                                                               + 1.0/dz_min/dz_min));
    int Rein_Step = 500;
    for (int Rein_step = 0; Rein_step <= Rein_Step; Rein_step++){
        //Fill LSBorder
        for (int lev = 0; lev <= finest_level; lev++){
            FillPatch(LSBorder[lev], lev, t_new[0],
                      LS, t_new,
                      LS, t_new,
                      0, 0, Num_Phase,
                      LS_bcs, 0,
                      &lincc_interp);
        }
        //Compute dLSdt
        for (int lev = 0; lev <= finest_level; lev++)
        {
            const auto dX = geom[lev].CellSizeArray();
            AMREX_D_TERM(Real dx = dX[0];,
                         Real dy = dX[1];,
                         Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                for (MFIter mfi(dLSdt[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx  = mfi.tilebox();
                
                    Array4<Real const> LS0Array   = LS0     [lev].const_array(mfi);
                    Array4<Real const> LSArray    = LSBorder[lev].const_array(mfi);
                    Array4<Real      > dLSdtArray = dLSdt   [lev].array(mfi);
                
                    amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        amrex::Real Diff = std::max({AMREX_D_DECL(dx, dy, dz)});
                        for (int i_Phase = 0; i_Phase < Num_Phase; ++i_Phase){
                            amrex::Real LS0_c = LS0Array(i,j,k,i_Phase);
                            amrex::Real Sign_LS0 = LS0_c/std::sqrt(LS0_c*LS0_c + Diff);

                            amrex::Real Laplace_LS;
                            FDM_Laplace_K(i, j, k, i_Phase,
                                          Laplace_LS,
                                          LSArray,
                             AMREX_D_DECL(dx, dy, dz));

                            amrex::Real Abs_Grad_LS = 0.0;
                            for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                                amrex::Real Grad_LS;
                                FDM_Gradient_CenterFromCenter_K(i, j, k, i_Phase,
                                                                Grad_LS,
                                                                LSArray,
                                                                dX[idim],
                                                                idim);
                                Abs_Grad_LS += Grad_LS*Grad_LS;
                            }
                            Abs_Grad_LS = std::sqrt(Abs_Grad_LS);

                            dLSdtArray(i,j,k,i_Phase) = -Sign_LS0*(Abs_Grad_LS - 1.0) + Diff*Laplace_LS;
                        } // end i_Phase
                    });
                } // end mfi
            } // end omp
        } // end lev
        //Update LS
        for (int lev = 0; lev <= finest_level; lev++){
            amrex::MultiFab::Saxpy(LS[lev], Rein_dt, dLSdt[lev], 0, 0, Num_Phase, 0);
        }
    } // end Rein_step
    //Fill LSBorder
    for (int lev = 0; lev <= finest_level; lev++){
        FillPatch(LSBorder[lev], lev, t_new[0],
                  LS, t_new,
                  LS, t_new,
                  0, 0, Num_Phase,
                  LS_bcs, 0,
                  &lincc_interp);
    }
}


// Compute Normal vector at cell centers from function f
void
Compressible_PhaseField::PhaseField_Normal(amrex::Vector<amrex::MultiFab>& NormalBorder,
                                           amrex::Vector<amrex::MultiFab> const& fBorder,
                                           int ID_Phase)
{
    // =======================================================
    // Compute Normal at all level
    // =======================================================
    //Define Normal
    Vector<MultiFab> Normal(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++){
        Normal[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }
    //Normal boundary conditions
    Vector<BCRec> Normal_bcs(AMREX_SPACEDIM);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        //Periodic
        if (geom[0].isPeriodic(idim)){
            for (int i_comp = 0; i_comp < AMREX_SPACEDIM; i_comp++){
                Normal_bcs[i_comp].setLo(idim, BCType::int_dir);
                Normal_bcs[i_comp].setHi(idim, BCType::int_dir);
            }
        }
        //Homogeneous Neumann
        else{
            for (int i_comp = 0; i_comp < AMREX_SPACEDIM; i_comp++){
                Normal_bcs[i_comp].setLo(idim, BCType::foextrap);
                Normal_bcs[i_comp].setHi(idim, BCType::foextrap);   
            }
        }
    }
    //Fill Normal
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const auto dX = geom[lev].CellSizeArray();
        AMREX_D_TERM(Real dx = dX[0];,
                     Real dy = dX[1];,
                     Real dz = dX[2]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(Normal[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx  = mfi.tilebox();
                
                Array4<Real const> fArray      = fBorder[lev].const_array(mfi);
                Array4<Real      > NormalArray = Normal [lev].array(mfi);
                
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    amrex::Array<amrex::Real,AMREX_SPACEDIM> Normal_f;
                    FDM_Normal_CenterFromCenter_K(i, j, k, ID_Phase,
                                                  Normal_f,
                                                  fArray,
                                     AMREX_D_DECL(dx, dy, dz));
                    for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                        NormalArray(i,j,k,idim) = Normal_f[idim];
                    }
                });
            } // end mfi
        } // end omp
    } // end lev
    //Fill NormalBorder
    for (int lev = 0; lev <= finest_level; lev++)
    {
        NormalBorder[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);
        FillPatch(NormalBorder[lev], lev, t_new[0],
                  Normal, t_new,
                  Normal, t_new,
                  0, 0, AMREX_SPACEDIM,
                  Normal_bcs, 0,
                  &lincc_interp);
    }
}


//Compute Phase-Field time step
void
Compressible_PhaseField::PhaseField_ComputeDt (Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, amrex::Real PF_time)
{
    Real PF_cfl = h_parm->PhaseField_Parm.cfl;

    Vector<Real> PF_dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        Real dx = std::min({AMREX_D_DECL(geom[lev].CellSize(0),geom[lev].CellSize(1),geom[lev].CellSize(2))});

        Real cmax = 0.0;
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
            for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
                cmax = std::max({ cmax, Grad_Q[lev][idim].norminf(i_Phase, 0, true) });
            } // end i_Phase
        } // end idim
        // cmax = PF_c_max[lev].norminf(0, 0, true);

        PF_dt_tmp[lev] = PF_cfl*dx/(cmax + 1.0e-14);
        if (PF_time == 0.0){
            PF_dt[lev] = PF_dt_tmp[lev];
        }
    } // end lev
    ParallelDescriptor::ReduceRealMin(&PF_dt_tmp[0], PF_dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real PF_dt_0 = PF_dt_tmp[0];
    int n_factor = 1;

    for (int lev = 0; lev <= finest_level; ++lev) {
        if (timestep_change_limiter && PF_dt[lev] > 0.0){
            PF_dt_tmp[lev] = std::min(PF_dt_tmp[lev], change_max*PF_dt[lev]);
        }
        n_factor *= nsubsteps[lev];
        PF_dt_0 = std::min(PF_dt_0, n_factor*PF_dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = (1.0e-3)*PF_dt_0;

    if (PF_time + PF_dt_0 > dt[0] - eps) {
        PF_dt_0 = dt[0] - PF_time;
    }

    PF_dt[0] = PF_dt_0;

    if (do_subcycle){
        for (int lev = 1; lev <= finest_level; ++lev) {
            PF_dt[lev] = PF_dt[lev-1] / nsubsteps[lev];
        }
    }
    else{
        amrex::Real PF_dt_min = PF_dt_0;
        for (int lev = 0; lev <= finest_level; ++lev) {
            PF_dt_min = std::min(PF_dt_min, PF_dt_tmp[lev]);
        }
        for (int lev = 0; lev <= finest_level; ++lev) {
            PF_dt[lev] = PF_dt_min;
        }
    }

}


//Explicit advancement in a time step
void 
Compressible_PhaseField::PhaseField_ExplicitTimeStep (amrex::Real cur_time,
                                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q)
{

    int lev = 0;
    int iteration = 1;
    int stage = 0;
    if (do_subcycle){
        PhaseField_TimeStepWithSubcycling(dof_new, 
                                          dof_old,
                                          Grad_Q,
                                          lev, 
                                          cur_time, 
                                          iteration, 
                                          stage);
    }
    else{
        PhaseField_TimeStepNoSubcycling(dof_new,
                                        dof_old,
                                        Grad_Q,
                                        cur_time,
                                        iteration);
    }
    
}



// Compute Grad_Q from Div_J
void
Compressible_PhaseField::PhaseField_GradQPF(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                                            amrex::Vector<amrex::MultiFab> const& Div_J,
                                            amrex::Vector<amrex::MultiFab> const& PFBorder,
                                            amrex::Real const& Eta,
                                            amrex::Real const& Eta_Factor,
                                            int ID_Weight,
                                            LinearSystem_Parameter const& LinearSystem_Parm)
{
    const int Num_Phase = PFBorder[0].nComp();

    // =======================================================
    // Define scalars
    // =======================================================
    amrex::Real ascalar =  0.0;
    amrex::Real bscalar = -1.0;
    
    // ==========================================================
    // Boundary condition
    // ==========================================================
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> LinOpBC_lo;
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> LinOpBC_hi;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        //Periodic
        if ( geom[0].isPeriodic(idim) ){
            LinOpBC_lo[idim] = LinOpBCType::Periodic;
            LinOpBC_hi[idim] = LinOpBCType::Periodic;
        }
        //Homogeneous Neumann
        else{
            LinOpBC_lo[idim] = LinOpBCType::Neumann;
            LinOpBC_hi[idim] = LinOpBCType::Neumann;
        }
    }
    int ID_CoarseFineBC = 0;
    amrex::MultiFab const* CoarseFineBCData = nullptr;
    int CoarseFineBCRatio = 2;
    
    // =======================================================
    // Define coefficients and solutions
    // =======================================================
    amrex::Vector<amrex::MultiFab> acoef(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> bcoef(finest_level+1);
    amrex::Vector<amrex::MultiFab> QPFi(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> Flux_QPFi(finest_level+1);
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> Grad_QPFi(finest_level+1);
    amrex::Vector<amrex::MultiFab> Div_Ji(finest_level+1);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        acoef[lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        acoef[lev].setVal(0.0);
        
        QPFi[lev] = MultiFab(grids[lev], dmap[lev], 1, 1);
        
        Div_Ji[lev] = MultiFab(grids[lev], dmap[lev], 1, 0);
        
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            bcoef[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
            Flux_QPFi[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
            Grad_QPFi[lev][idim] = MultiFab(ba, dmap[lev], 1, 0);
        }
    }

    // =======================================================
    // Compute Grad_QPF for all phases
    // =======================================================
    for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
        // =======================================================
        // Compute bcoef at all level
        // =======================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                amrex::Real ep = 1.0e-3;
                for (MFIter mfi(Div_J[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();

                    Array4<Real const> PFArray = PFBorder[lev].const_array(mfi);
                    
                    for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                        Array4<Real> bcoefArray = bcoef[lev][idim].array(mfi);

                        // the faces of this tile only, so that no two OpenMP threads write one face
                        amrex::ignore_unused(bx);
                        const Box bx_face = mfi.nodaltilebox(idim);
                        amrex::ParallelFor(bx_face,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k)
                        {
                            amrex::Real PF_c = PFArray(i,j,k,i_Phase);
                            amrex::Real PF_m;
                            //x
                            if (idim == 0){
                                PF_m = PFArray(i-1,j,k,i_Phase);
                            }
                            //y
                            else if (idim == 1){
                                PF_m = PFArray(i,j-1,k,i_Phase);
                            }
                            //z
                            else if (idim == 2){
                                PF_m = PFArray(i,j,k-1,i_Phase);
                            }
                            
                            int ID_Derivative = 0;
                            amrex::Real Psi_c;
                            Mapping_Backward(Psi_c,
                                             PF_c,
                                             Eta,
                                             ID_Derivative);
                            amrex::Real Psi_m;
                            Mapping_Backward(Psi_m,
                                             PF_m,
                                             Eta,
                                             ID_Derivative);
                            Mapping_Forward(PF_c,
                                            Psi_c,
                                            Eta_Factor*Eta,
                                            ID_Derivative);
                            Mapping_Forward(PF_m,
                                            Psi_m,
                                            Eta_Factor*Eta,
                                            ID_Derivative);
                            amrex::Real PF_Face;
                            Interpolation2Face_Linear_D(PF_Face,
                                                        PF_m,
                                                        PF_c);
                            amrex::Real Weight_PF;
                            WeightFunction_PhaseField(Weight_PF,
                                                      PF_Face,
                                                      ID_Derivative, ID_Weight);
                            bcoefArray(i,j,k,0) = Weight_PF + ep;
                        });
                    } // end idim
                } // end mfi
            } // end omp
        } // end lev
        // ==========================================================
        // Average down bcoef
        // ==========================================================
        AverageDownFaces(bcoef);
        
        // ==========================================================
        // Compute Grad_QPFi
        // ==========================================================
        for (int lev = 0; lev <= finest_level; lev++)
        {
            QPFi[lev].setVal(0.0);
            MultiFab::Copy(Div_Ji[lev], Div_J[lev], i_Phase, 0, 1, 0);
        }
        SolveABecLaplacianAllLevels_HomBC(QPFi, Flux_QPFi, Grad_QPFi,
                                          Geom(0, finest_level), grids, dmap, ref_ratio,
                                          ascalar, bscalar, acoef, bcoef, Div_Ji,
                                          LinOpBC_lo, LinOpBC_hi,
                                          ID_CoarseFineBC, CoarseFineBCData, CoarseFineBCRatio,
                                          LinearSystem_Parm);
        // ===============================================================
        // Average down Grad_QPFi
        // ===============================================================
        AverageDownFaces(Grad_QPFi);
        
        // ===============================================================
        // Copy Grad_QPFi to Grad_QPF
        // ===============================================================
        for (int lev = 0; lev <= finest_level; lev++){
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                MultiFab::Copy(Grad_QPF[lev][idim], Grad_QPFi[lev][idim], 0, i_Phase, 1, 0);
            } // end idim
        } // end lev
        
    } // end i_Phase
    
}


// Compute Grad_QPF from J_Flux
void
Compressible_PhaseField::PhaseField_GradQPF(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                                            amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> const& J_Flux,
                                            amrex::Vector<amrex::MultiFab> const& PFBorder,
                                            amrex::Real const& Eta,
                                            amrex::Real const& Eta_Factor)
{
    const int Num_Phase = PFBorder[0].nComp();
    // =======================================================
    // Fill Grad_QPF for all phases
    // =======================================================
    amrex::Real epsilon = 1.0e-3;
    for (int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            for (MFIter mfi(PFBorder[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                Array4<Real const> PFArray = PFBorder[lev].const_array(mfi);
                
                for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
                {
                    Array4<Real const> J_FluxArray   = J_Flux  [lev][idim].const_array(mfi);
                    Array4<Real      > Grad_QPFArray = Grad_QPF[lev][idim].array(mfi);
                    
                    // the faces of this tile only, so that no two OpenMP threads write one face
                    amrex::ignore_unused(bx);
                    const Box bx_face = mfi.nodaltilebox(idim);
                    amrex::ParallelFor(bx_face,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        for (int i_Phase = 0; i_Phase < Num_Phase; i_Phase++){
                            amrex::Real PF_c = PFArray(i,j,k,i_Phase);
                            amrex::Real PF_m;
                            //x
                            if (idim == 0){
                                PF_m = PFArray(i-1,j,k,i_Phase);
                            }
                            //y
                            else if (idim == 1){
                                PF_m = PFArray(i,j-1,k,i_Phase);
                            }
                            //z
                            else if (idim == 2){
                                PF_m = PFArray(i,j,k-1,i_Phase);
                            }
                            
                            int ID_Derivative = 0;
                            amrex::Real Psi_c;
                            Mapping_Backward(Psi_c,
                                             PF_c,
                                             Eta,
                                             ID_Derivative);
                            amrex::Real Psi_m;
                            Mapping_Backward(Psi_m,
                                             PF_m,
                                             Eta,
                                             ID_Derivative);
                            
                            Mapping_Forward(PF_c,
                                            Psi_c,
                                            Eta_Factor*Eta,
                                            ID_Derivative);
                            Mapping_Forward(PF_m,
                                            Psi_m,
                                            Eta_Factor*Eta,
                                            ID_Derivative);
                            
                            amrex::Real PF_Face;
                            Interpolation2Face_Linear_D(PF_Face,
                                                        PF_m,
                                                        PF_c);
                            
                            amrex::Real PF_Flux = J_FluxArray(i,j,k,i_Phase);
                            amrex::Real GradQPF;
                            PhaseFieldFlux2GradQPF(GradQPF,
                                                   PF_Flux, PF_Face,
                                                   epsilon);
                            Grad_QPFArray(i,j,k,i_Phase) = GradQPF;
                        } // end i_Phase
                    });
                } // end idim
            } // end mfi
        } // end omp
    } // end lev
}


// Compute Grad_QPF from J_Flux with degenerate mobility
void
Compressible_PhaseField::PhaseField_GradQPF(amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                                            amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> const& J_Flux)
{
    // =======================================================
    // Fill Grad_QPF for all phases
    // =======================================================
    for (int lev = 0; lev <= finest_level; lev++)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
        {
            amrex::MultiFab::Copy(Grad_QPF[lev][idim], J_Flux[lev][idim], 0, 0, J_Flux[lev][idim].nComp(), 0);
        } // end idim
    } // end lev
}


//Euler advancement
void 
Compressible_PhaseField::PhaseField_ForwardEuler (Vector<MultiFab>& mf_new, 
                                                  Vector<MultiFab>& mf_old, 
                                                  Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                                  int lev, 
                                                  Real time, 
                                                  Real dt_lev, 
                                                  int ncycle)
{
    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    const int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }


#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }

    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif

    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);

    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_PreTimeStep(Uborder,ncycle,time,lev);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev],
                               lev, 
                               dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time,lev);
    PhaseField_PostTimeStep(U_new,ncycle,time,lev);

}

//RK2 advancement
void 
Compressible_PhaseField::PhaseField_SecondOrderSSPRK (Vector<MultiFab>& mf_new, 
                                                      Vector<MultiFab>& mf_old, 
                                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                                      int lev, 
                                                      Real time, 
                                                      Real dt_lev, 
                                                      int ncycle)
{
    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    const int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }

#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }
    
    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif

    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);


    // RK2 stage 1
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_PreTimeStep(Uborder,ncycle,time,lev);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               Real(0.5)*dt_lev, 
                               ncycle, 
                               stage,
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time+dt_lev,lev);


    /* RK2 stage 2 */
    // After fillpatch Uborder = U^n+dt*dUdt^n = U^*
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time+dt_lev,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               Real(0.5)*dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    // U_new = 0.5*(Uborder+U_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(U_new, Real(0.5), Uborder, 0, Real(0.5), U_old, 0, 0, NSTATE, 0);
    // U_new += 0.5*dt*dUdt
    MultiFab::Saxpy(U_new, Real(0.5)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*

    PhaseField_PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PhaseField_PostTimeStep(U_new,ncycle,time+dt_lev,lev);

}

//RK3 advancement
void 
Compressible_PhaseField::PhaseField_ThirdOrderSSPRK (Vector<MultiFab>& mf_new, 
                                                     Vector<MultiFab>& mf_old, 
                                                     Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                                     int lev, 
                                                     Real time, 
                                                     Real dt_lev, 
                                                     int ncycle)
{

    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    const int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }


#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }

    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif


    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);


    // RK3 stage 1
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_PreTimeStep(Uborder,ncycle,time,lev);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev],
                               lev,
                               (Real(1.0)/Real(6.0))*dt_lev,
                               ncycle,
                               stage,
                               fr_as_crse,
                               fr_as_fine,
                               fr_as_crse_nc,
                               fr_as_fine_nc);
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time+dt_lev,lev);


    /* RK3 stage 2 */
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time+dt_lev,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev],
                               lev,
                               (Real(1.0)/Real(6.0))*dt_lev,
                               ncycle,
                               stage,
                               fr_as_crse,
                               fr_as_fine,
                               fr_as_crse_nc,
                               fr_as_fine_nc);
    MultiFab::LinComb(U_new, Real(0.25), Uborder, 0, Real(0.75), U_old, 0, 0, NSTATE, 0);
    MultiFab::Saxpy(U_new, Real(0.25)*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time+Real(0.5)*dt_lev,lev);


    /* RK3 stage 3 */
    // U_new = U^(2) is at time+dt/2, not at t_new[lev] = time+dt, so fill it as a stage state
    stage = 2;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_new,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev],
                               lev,
                               (Real(2.0)/Real(3.0))*dt_lev,
                               ncycle,
                               stage,
                               fr_as_crse,
                               fr_as_fine,
                               fr_as_crse_nc,
                               fr_as_fine_nc);
    MultiFab::LinComb(U_new, (Real(2.0)/Real(3.0)), Uborder, 0, (Real(1.0)/Real(3.0)), U_old, 0, 0, NSTATE, 0);
    MultiFab::Saxpy(U_new, (Real(2.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PhaseField_PostTimeStep(U_new,ncycle,time+dt_lev,lev);

}

//RK4 advancement
void 
Compressible_PhaseField::PhaseField_FourthOrderRK (Vector<MultiFab>& mf_new, 
                                                   Vector<MultiFab>& mf_old, 
                                                   Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                                   int lev, 
                                                   Real time, 
                                                   Real dt_lev, 
                                                   int ncycle)
{

    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    const int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }


#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }

    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif


    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab U_temp     (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);

    amrex::MultiFab::Copy(U_temp, U_new,0,0,NSTATE,0);


    //=====================================================================
    //RK4 stage 1 
    //=====================================================================
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PhaseField_PreTimeStep(Uborder,ncycle,time,lev);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               (Real(1.0)/Real(6.0))*dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    /* u_1 = u_n + 1/2 dt R(u_n) */
    MultiFab::LinComb(U_temp, Real(1.0), Uborder, 0, Real(0.5)*dt_lev,            dUdt, 0, 0, NSTATE, 0);
    /* u_* = u_n + 1/6 dt R(u_n) */
    MultiFab::LinComb(U_new,  Real(1.0), Uborder, 0, (Real(1.0)/Real(6.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_temp,ncycle,time+Real(0.5)*dt_lev,lev);


    //=====================================================================
    //RK4 stage 2 
    //=====================================================================
    // the stage state u_1 at time+dt/2 is in U_temp, while U_new holds the running sum u_*
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               (Real(1.0)/Real(3.0))*dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    /* u_2 = u_n + 1/2 dt R(u_1) */
    MultiFab::LinComb(U_temp, Real(1.0), U_old, 0, Real(0.5)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    /* u_** = u_* + 1/3 dt R(u_1) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_temp,ncycle,time+Real(0.5)*dt_lev,lev);


    //=====================================================================
    //RK4 stage 3 
    //=====================================================================
    // the stage state u_2 at time+dt/2 is in U_temp
    stage = 2;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               (Real(1.0)/Real(3.0))*dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    /* u_3 = u_n + dt R(u_2) */
    MultiFab::LinComb(U_temp, Real(1.0), U_old, 0, Real(1.0)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    /* u_*** = u_** + 1/3 dt R(u_2) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_temp,ncycle,time+dt_lev,lev);


    //=====================================================================
    //RK4 stage 4
    //=====================================================================
    // the stage state u_3 at time+dt is in U_temp, not in U_new, so it is filled as a stage state
    stage = 3;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    PhaseField_compute_dUdt_FV(dUdt,
                               Uborder,
                               Grad_Q[lev], 
                               lev, 
                               (Real(1.0)/Real(6.0))*dt_lev, 
                               ncycle, 
                               stage, 
                               fr_as_crse, 
                               fr_as_fine, 
                               fr_as_crse_nc, 
                               fr_as_fine_nc);
    /* u_n+1 = u_*** + 1/6 dt R(u_3) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(6.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    PhaseField_PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PhaseField_PostTimeStep(U_new,ncycle,time+dt_lev,lev);

}



// Advance a level by PF_dt (includes a recursive call for finer levels)
void
Compressible_PhaseField::PhaseField_TimeStepWithSubcycling (Vector<MultiFab>& mf_new,
                                                            Vector<MultiFab>& mf_old,
                                                            Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                                            int lev, 
                                                            amrex::Real time, 
                                                            int iteration,
                                                            int stage)
{

    BL_PROFILE("PhaseField_timeStepWithSubcycling()");

    // if (stage == 0){
    //     if (Verbose()) {
    //         // amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
    //         amrex::Print() << "[PF] Phase Field ADVANCE with time = " << t_new[lev]
    //                        << " PF_dt = " << PF_dt[lev] << std::endl;
    //     }
    // }


    t_old[lev] = t_new[lev];
    t_new[lev] += PF_dt[lev];

    // Advance a single level for a single time step, and update flux registers
    PhaseField_AdvanceAtLevel(mf_new, 
                              mf_old, 
                              Grad_Q, 
                              lev, 
                              time, 
                              PF_dt[lev], 
                              nsubsteps[lev]);



    // if (stage == 0){
    //     // ++isubstep;

    //     if (Verbose())
    //     {
    //         amrex::Print() << "[Level " << lev << " step " << isubstep[lev] << "] ";
    //         amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    //     }
    // }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            PhaseField_TimeStepWithSubcycling(mf_new, 
                                              mf_old, 
                                              Grad_Q,
                                              lev+1, 
                                              time+(i-1)*PF_dt[lev+1], 
                                              i, 
                                              stage);
        }

        if (do_reflux)
        {
            PhaseField_RefluxLev(lev);
        }

        AverageDownTo(mf_new,lev); // average lev+1 down to lev
    }

}


// Advance all the levels with the same PF_dt
void
Compressible_PhaseField::PhaseField_TimeStepNoSubcycling (Vector<MultiFab>& mf_new,
                                                          Vector<MultiFab>& mf_old,
                                                          Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                                          Real time,
                                                          int iteration)
{
    PhaseField_AdvanceAllLevels(mf_new,
                                mf_old,
                                Grad_Q,
                                time,
                                PF_dt[0],
                                iteration);

    // Make sure the coarser levels are consistent with the finer levels
    AverageDown (mf_new);

    // istep counts the time steps of the flow (timeStepNoSubcycling); the Phase-Field steps inside
    // one of them do not count, as in PhaseField_TimeStepWithSubcycling
}


void
Compressible_PhaseField::PhaseField_RefluxLev (int lev)
{
    // update lev based on coarse-fine flux mismatch
    flux_reg[lev+1]->Reflux(dof_new[lev], 1.0, 0, 0, NSTATE, geom[lev]);
}


void 
Compressible_PhaseField::SetPhysicsBC ()
{
    bcs.resize(NSTATE);

    //
    // Components are:
    //  Interior,        Inflow,          Outflow,          Symmetry,              SlipWall,            NoSlipWall
    //
    static int scalar_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
    };

    static int norm_vel_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_odd,  BCType::reflect_odd,  BCType::reflect_odd
    };

    static int tang_vel_bc[] =
    {
        BCType::int_dir, BCType::ext_dir, BCType::foextrap, BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd
    };




    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
/*
        bcs[INDEX_Mass1].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Mass2].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_Energy].setLo(idim, scalar_bc[lo_bc[idim]]);
        bcs[INDEX_VolumeFraction1].setLo(idim, scalar_bc[lo_bc[idim]]);

        bcs[INDEX_Mass1].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_Mass2].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_Energy].setHi(idim, scalar_bc[hi_bc[idim]]);
        bcs[INDEX_VolumeFraction1].setHi(idim, scalar_bc[hi_bc[idim]]);
*/
    }


#if (AMREX_SPACEDIM == 1)
    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);

#elif (AMREX_SPACEDIM == 2)
    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumX].setLo(1, tang_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumX].setHi(1, tang_vel_bc[hi_bc[1]]);

    bcs[INDEX_MomentumY].setLo(0, tang_vel_bc[lo_bc[0]]);    
    bcs[INDEX_MomentumY].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumY].setLo(1, norm_vel_bc[lo_bc[1]]);    
    bcs[INDEX_MomentumY].setHi(1, norm_vel_bc[hi_bc[1]]);

#elif (AMREX_SPACEDIM == 3)

    bcs[INDEX_MomentumX].setLo(0, norm_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumX].setHi(0, norm_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumX].setLo(1, tang_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumX].setHi(1, tang_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumX].setLo(2, tang_vel_bc[lo_bc[2]]);
    bcs[INDEX_MomentumX].setHi(2, tang_vel_bc[hi_bc[2]]);

    bcs[INDEX_MomentumY].setLo(0, tang_vel_bc[lo_bc[0]]);
    bcs[INDEX_MomentumY].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumY].setLo(1, norm_vel_bc[lo_bc[1]]);
    bcs[INDEX_MomentumY].setHi(1, norm_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumY].setLo(2, tang_vel_bc[lo_bc[2]]);
    bcs[INDEX_MomentumY].setHi(2, tang_vel_bc[hi_bc[2]]);

    bcs[INDEX_MomentumZ].setLo(0, tang_vel_bc[lo_bc[0]]);    
    bcs[INDEX_MomentumZ].setHi(0, tang_vel_bc[hi_bc[0]]);
    bcs[INDEX_MomentumZ].setLo(1, tang_vel_bc[lo_bc[1]]);    
    bcs[INDEX_MomentumZ].setHi(1, tang_vel_bc[hi_bc[1]]);
    bcs[INDEX_MomentumZ].setLo(2, norm_vel_bc[lo_bc[2]]);    
    bcs[INDEX_MomentumZ].setHi(2, norm_vel_bc[hi_bc[2]]);
#endif

}


// Face MultiFabs of level lev: the NSTATE conservative fluxes and the NC_TERMS face quantities of
// the non-conservative terms (not defined without NONCONSERVATIVE)
void
Compressible_PhaseField::DefineFaceFluxes (Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                           Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                           int lev)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        BoxArray ba = grids[lev];
        ba.surroundingNodes(idim);
        fluxes[idim].define(ba, dmap[lev], NSTATE, 0);
#if (NONCONSERVATIVE == true)
        fluxes_nc[idim].define(ba, dmap[lev], NC_TERMS, 0);
#else
        amrex::ignore_unused(fluxes_nc);
#endif
    }
}


// Face fluxes of level lev from the state U with ghost cells: in fluxes the conservative fluxes per
// unit area (hyperbolic minus diffusive), in fluxes_nc the face quantities of the
// non-conservative terms, and the wave speeds in c_max[lev] (reset at stage 0). time is the time of
// the stage state. Each tile computes the faces of nodaltilebox, which no other tile has, and the
// c_max of a cell is written only from its own faces (i,j,k), which belong to the same tile: no
// two OpenMP threads write the same face or cell
void
Compressible_PhaseField::ComputeFaceFluxes (MultiFab const& U,
                                            Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                            Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                            int lev,
                                            Real time,
                                            int stage)
{
BL_PROFILE("ComputeFaceFluxes()");

    MultiFab& c_max_lev = c_max[lev];
    if (stage == 0){
        c_max_lev.setVal(0.0);
    }

    auto const prob_lo = Geom(lev).ProbLoArray();
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
    amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
    amrex::ignore_unused(prob_lo, time, dy, dz);

    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(U,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Pull the data into an array
            Array4<Real const> statein   = U.const_array(mfi);
            Array4<Real      > c_max_new = c_max_lev.array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                Array4<Real> flux = fluxes[idim].array(mfi);
#if (NONCONSERVATIVE == true)
                Array4<Real> fluxNC = fluxes_nc[idim].array(mfi);
#endif
{BL_PROFILE("ComputeFaceFluxes::{ idim-direction ParallelFor }");
                // the faces of this tile only (the high face of a box with its last tile)
                const Box bx_face = mfi.nodaltilebox(idim);
                amrex::ParallelFor(bx_face,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Array<Real,NSTATE> fhat;
                    Array<Real,NC_TERMS> fhatNC;
                    Real cmax;
                    FVM_Conservative2FluxHyperbolic_K(i, j, k,
                                                      fhat, fhatNC, cmax,
                                                      statein, *lparm,
                                                      dX,
                                                      idim,
                                                      lparm->FiniteVolume_Parm, lparm->Physics_Parm);
                    c_max_new(i,j,k,0) = std::max(cmax, c_max_new(i,j,k,0));
#if (DIFFUSION == true)
                    Array<Real,NSTATE> fhatD;
                    amrex::Real dmax;
#if (PHYSICS == SIXEQS)
                    // SIXEQS also writes the viscous energy flux to a non-conservative component
                    FDM_Conservative2FluxDiffusion_K(i, j, k,
                                                     fhatD, fhatNC, dmax,
                                                     statein,
                                                     prob_lo,
                                        AMREX_D_DECL(dx,dy,dz),
                                                     time,
                                                     *lparm,
                                                     idim);
#else
                    FDM_Conservative2FluxDiffusion_K(i, j, k,
                                                     fhatD, dmax,
                                                     statein,
                                                     prob_lo,
                                        AMREX_D_DECL(dx,dy,dz),
                                                     time,
                                                     *lparm,
                                                     idim);
#endif
                    c_max_new(i,j,k,1) = std::max(dmax, c_max_new(i,j,k,1));
#endif
                    for (int iState = 0; iState < NSTATE; iState++){
                        flux(i,j,k,iState) = fhat[iState];
#if (DIFFUSION == true)
                        flux(i,j,k,iState) += -fhatD[iState];
#endif
                    }
#if (NONCONSERVATIVE == true)
                    for (int iState = 0; iState < NC_TERMS; iState++){
                        fluxNC(i,j,k,iState) = fhatNC[iState];
                    }
#endif
                });
}
            } // end idim
        } // end mfi
    } // end omp
}


// Right-hand side dUdt of level lev (overwritten) from the face fluxes of ComputeFaceFluxes, possibly
// with the averaged fluxes of the finer level on the covered faces: the conservative surface
// integral, the user source term and the non-conservative terms, which read the same face
// MultiFabs. U is the state with ghost cells the fluxes were computed from, time its time; dt_lev
// is passed to the user source term
void
Compressible_PhaseField::FluxDivergence (MultiFab& dUdt_mf,
                                         MultiFab const& U,
                                         Array<MultiFab,AMREX_SPACEDIM> const& fluxes,
                                         Array<MultiFab,AMREX_SPACEDIM> const& fluxes_nc,
                                         int lev,
                                         Real time,
                                         Real dt_lev)
{
BL_PROFILE("FluxDivergence()");

    auto const prob_lo = Geom(lev).ProbLoArray();
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
    amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
    amrex::ignore_unused(time, dt_lev, dy, dz, fluxes_nc);

    Parm const* lparm = d_parm;

    const int lID_COORDSYS = ID_COORDSYS;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(dUdt_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            Array4<Real const> statein = U.const_array(mfi);
            Array4<Real      > dUdt    = dUdt_mf.array(mfi);

            AMREX_D_TERM(Array4<Real const> fluxx_c = fluxes[0].const_array(mfi);,
                         Array4<Real const> fluxy_c = fluxes[1].const_array(mfi);,
                         Array4<Real const> fluxz_c = fluxes[2].const_array(mfi));
{BL_PROFILE("COMPAS::compute_dUdt_FV::surface_integral()");
            // Do a conservative update
            // Forward Euler
            // ===========================
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (lID_COORDSYS == 0){
                    for (int iState = 0; iState < NSTATE; iState++){
                        FVM_SurfaceIntegral(i, j, k, iState,
                                         dUdt, statein,
                            AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
                            AMREX_D_DECL(dx,dy,dz));
                    }
                }
                else {
                    FVM_SurfaceIntegral_Coord(i,j,k,
                                           dUdt,statein,
                              AMREX_D_DECL(fluxx_c,fluxy_c,fluxz_c),
                                           lID_COORDSYS,
                                           dX,
                                           prob_lo);
                }
            });
}

#ifdef USER_SOURCE_TERM
{BL_PROFILE("COMPAS::compute_dUdt_FV::user_source_term()");
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real x = prob_lo[0] + (Real(i)+0.5)*dx;
#if (AMREX_SPACEDIM > 1)
                amrex::Real y = prob_lo[1] + (Real(j)+0.5)*dy;
#else
                amrex::Real y = 0.0;
#endif
#if (AMREX_SPACEDIM > 2)
                amrex::Real z = prob_lo[2] + (Real(k)+0.5)*dz;
#else
                amrex::Real z = 0.0;
#endif
                amrex::Real t = time;

                amrex::Array<amrex::Real, NSTATE> Residual;
                for (int iState = 0; iState < NSTATE; iState++){
                    Residual[iState] = 0.0;
                }

                user_source_term(
                    Residual,
                    statein,
                    i, j, k,
                    x, y, z, t, dt_lev,
                    AMREX_D_DECL(dx,dy,dz),
                    *lparm);

                for (int iState = 0; iState < NSTATE; iState++){
                    dUdt(i,j,k,iState) += Residual[iState];
                }
            });
}
#endif

#if (NONCONSERVATIVE == true)
            AMREX_D_TERM(Array4<Real const> fluxxNC_c = fluxes_nc[0].const_array(mfi);,
                         Array4<Real const> fluxyNC_c = fluxes_nc[1].const_array(mfi);,
                         Array4<Real const> fluxzNC_c = fluxes_nc[2].const_array(mfi));
{BL_PROFILE("COMPAS::compute_dUdt_FV::surface_integral_nc()");
            // Add the non-conservative contribution
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FVM_SurfaceIntegral_NC(i, j, k,
                                    dUdt, statein,
                       AMREX_D_DECL(fluxx_c  , fluxy_c  , fluxz_c  ),
                       AMREX_D_DECL(fluxxNC_c, fluxyNC_c, fluxzNC_c),
                       AMREX_D_DECL(dx       , dy       , dz       ),
                                    *lparm,
                                    lparm->Physics_Parm);
            });
}
#endif
        } // end mfi
    } // end omp
}


// Abort if the wave speeds of level lev give c_max dt_lev > dx in a direction (run.check_cfl)
void
Compressible_PhaseField::CheckCFL (int lev, Real dt_lev)
{
    // ======== CFL CHECK, MOVED OUTSIDE MFITER LOOP =========
/* TODO: Make sure it is ok to not do this on the GPU/see if is actually that expensive*/
    if (check_cfl){
#if (!AMREX_USE_GPU)
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
        Real cmax = c_max[lev].norminf(0,0,true);

        if (AMREX_D_TERM(cmax*dt_lev > dx, ||
                         cmax*dt_lev > dy, ||
                         cmax*dt_lev > dz))
        {
            amrex::AllPrint() << "cmax = " << cmax
                              << ", dt = " << dt_lev << " dx = " << dx << " " << dy << " " << dz << std::endl;
            amrex::Abort("CFL violation. use smaller adv.cfl.");
        }
#else
        amrex::ignore_unused(lev, dt_lev);
#endif
    }
}


// Multiply the face fluxes of level lev by the face areas, for the flux registers. On whole face
// MultiFabs after the face and update loops: each face is multiplied once, and no tile reads a
// face that another tile scales
void
Compressible_PhaseField::ScaleFaceFluxesByArea (Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                                Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                                int lev)
{
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
    amrex::ignore_unused(fluxes_nc);
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
        amrex::Real dArea;
        if (idim == 0){
            dArea = dy*dz;
        }
        else if (idim == 1){
            dArea = dx*dz;
        }
        else {
            dArea = dx*dy;
        }
        fluxes[idim].mult(dArea, 0, NSTATE, 0);
#if (NONCONSERVATIVE == true)
        fluxes_nc[idim].mult(dArea, 0, NC_TERMS, 0);
#endif
    }
}


// Right-hand side of one level with its own face fluxes, for the per-level integrators of the
// subcycled time step; with do_reflux the area-weighted fluxes go to the flux registers
// (dt_lev = stage weight times dt of the level)
void
Compressible_PhaseField::compute_dUdt_FV (MultiFab& mf_new,
                                          MultiFab& mf_old,
                                          int lev,
                                          Real time,
                                          Real dt_lev,
                                          int ncycle,
                                          int stage,
                                          FluxRegister* fr_as_crse,
                                          FluxRegister* fr_as_fine,
                                          FluxRegister* fr_as_crse_nc,
                                          FluxRegister* fr_as_fine_nc)
{
BL_PROFILE("compute_dUdt_FV()");
    amrex::ignore_unused(ncycle, fr_as_crse_nc, fr_as_fine_nc);

    // construct the fluxes
    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    Array<MultiFab,AMREX_SPACEDIM> fluxes_nc;
    DefineFaceFluxes(fluxes, fluxes_nc, lev);

    // time is the time of the stage state, not t_old[lev], so that a time-dependent
    // source term is evaluated at the stage time
    ComputeFaceFluxes(mf_old, fluxes, fluxes_nc, lev, time, stage);

    FluxDivergence(mf_new, mf_old, fluxes, fluxes_nc, lev, time, dt_lev);

    CheckCFL(lev, dt_lev);

    if (do_reflux)
    {
        ScaleFaceFluxesByArea(fluxes, fluxes_nc, lev);
    }

    // increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // COMPAS solves U_t = -div(+F) for the conserved state U
    // The fluxes contain, e.g., F_{i+1/2,j} = F(U)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (fr_as_crse) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = -dt_lev;
            fr_as_crse->CrseInit(fluxes[idim], idim, 0, 0, NSTATE, scale, FluxRegister::ADD);
#if (NONCONSERVATIVE == true)
            // fr_as_crse_nc->CrseInit(fluxes_nc[idim], idim, 0, 0, NC_FLUX_COMP, scale, FluxRegister::ADD);
            fr_as_crse_nc->CrseInit(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale, FluxRegister::ADD);
#endif
        }
    }

    if (fr_as_fine) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = dt_lev;
            fr_as_fine->FineAdd(fluxes[idim], idim, 0, 0, NSTATE, scale);
#if (NONCONSERVATIVE == true)
            // fr_as_fine_nc->FineAdd(fluxes_nc[idim], idim, 0, 0, NC_FLUX_COMP, scale);
            fr_as_fine_nc->FineAdd(fluxes_nc[idim], idim, 0, 0, NC_TERMS, scale);
#endif
        }
    }
}


void 
Compressible_PhaseField::MaxVelocity (amrex::Real & max_V,
                                      amrex::MultiFab const& Conservative)
{

    max_V = 0.0;
    const BoxArray& ba = Conservative.boxArray();
    const DistributionMapping& dm = Conservative.DistributionMap();
    MultiFab Velocity(ba, dm, 1, 0);
    /* TODO: run on GPU*/
    Velocity.setVal(0.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(Conservative,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real const> ConservativeArray = Conservative.const_array(mfi);
            Array4<Real>       VelocityArray = Velocity.array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Array<Real,NSTATE> Conservative_c;
                for (int i_State = 0; i_State < NSTATE; i_State++){
                    Conservative_c[i_State] = ConservativeArray(i,j,k,i_State);
                }
                AMREX_D_TERM(Real u;, Real v;, Real w);
                Conservative2Velocity(AMREX_D_DECL(u, v, w),
                                                   Conservative_c);
                Real V = std::sqrt( AMREX_D_TERM(u*u, + v*v, + w*w) );
                VelocityArray(i,j,k) = V;
            });
        } // end mfi
    } // end omp
    max_V = std::max(max_V, Velocity.norminf(0,0,false));
}


void 
Compressible_PhaseField::MaxVelocity (amrex::Real & max_V,
                                      amrex::Vector<amrex::MultiFab> const& Conservative)
{
    max_V = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Real max_V_lev;
        MaxVelocity(max_V_lev,
                    Conservative[lev]);
        max_V = std::max(max_V, max_V_lev);
    } // end lev
}


void 
Compressible_PhaseField::MinEnergyBound (amrex::Real & min_EnergyBound,
                                         amrex::MultiFab const& Conservative,
                                         Parm const& EOS_Parameter)
{
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    min_EnergyBound = inf;
    const BoxArray& ba = Conservative.boxArray();
    const DistributionMapping& dm = Conservative.DistributionMap();
    MultiFab EnergyBound(ba, dm, 1, 0);
    /* TODO: run on GPU*/
    EnergyBound.setVal(inf);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Conservative,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real const> ConservativeArray = Conservative.const_array(mfi);
            Array4<Real>       EnergyBoundArray  = EnergyBound.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Array<Real,NSTATE> Conservative_c;
                for (int i_State = 0; i_State < NSTATE; i_State++){
                    Conservative_c[i_State] = ConservativeArray(i,j,k,i_State);
                }
                amrex::Real EB;
                Conservative2EnergyBound(EB,
                                         Conservative_c,
                                         EOS_Parameter);
                EnergyBoundArray(i,j,k,0) = EB;
            });
        } // end mfi
    } // end omp
    min_EnergyBound = std::min(min_EnergyBound, EnergyBound.min(0));
}


void 
Compressible_PhaseField::MinEnergyBound (amrex::Real & min_EnergyBound,
                                         amrex::Vector<amrex::MultiFab> const& Conservative,
                                         Parm const& EOS_Parameter)
{
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    min_EnergyBound = inf;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Real min_EnergyBound_lev;
        MinEnergyBound(min_EnergyBound_lev,
                       Conservative[lev],
                       EOS_Parameter);
        min_EnergyBound = std::min(min_EnergyBound, min_EnergyBound_lev);
    } // end lev
}


void 
Compressible_PhaseField::MinPressure (amrex::Real & min_Pressure,
                                      amrex::MultiFab const& Conservative,
                                      Parm const& EOS_Parameter)
{
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    min_Pressure = inf;
    const BoxArray& ba = Conservative.boxArray();
    const DistributionMapping& dm = Conservative.DistributionMap();
    MultiFab Pressure(ba, dm, 1, 0);
    /* TODO: run on GPU*/
    Pressure.setVal(inf);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Conservative,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real const> ConservativeArray = Conservative.const_array(mfi);
            Array4<Real>       PressureArray     = Pressure.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Array<Real,NSTATE> Conservative_c;
                for (int i_State = 0; i_State < NSTATE; i_State++){
                    Conservative_c[i_State] = ConservativeArray(i,j,k,i_State);
                }
                amrex::Real p;
                Conservative2Pressure(p,
                                      Conservative_c,
                                      EOS_Parameter);
                PressureArray(i,j,k,0) = p;
            });
        } // end mfi
    } // end omp
    min_Pressure = std::min(min_Pressure, Pressure.min(0));
}


void 
Compressible_PhaseField::MinPressure (amrex::Real & min_Pressure,
                                      amrex::Vector<amrex::MultiFab> const& Conservative,
                                      Parm const& EOS_Parameter)
{
    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    min_Pressure = inf;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Real min_Pressure_lev;
        MinPressure(min_Pressure_lev,
                    Conservative[lev],
                    EOS_Parameter);
        min_Pressure = std::min(min_Pressure, min_Pressure_lev);
    } // end lev
}

void 
Compressible_PhaseField::ForwardEuler (Vector<MultiFab>& mf_new, 
                                       Vector<MultiFab>& mf_old, 
                                       int lev, 
                                       Real time, 
                                       Real dt_lev, 
                                       int ncycle)
{
    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }


#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }
    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif

    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);
    MultiFab flagfab    (grids[lev],dmap[lev],1,1);
    Real min_EnergyBound;

    int stage = 0;
    dUdt.setVal(Real(0.0));

    FillPatch(Uborder, lev, time,
               mf_new, t_new,
               mf_old, t_old,
               0, 0, NSTATE);

    PreTimeStep(Uborder,ncycle,time,lev);

    compute_dUdt_FV(dUdt,Uborder,lev, time, dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_new,ncycle,time,lev);
    PostTimeStep(U_new,ncycle,time,lev);

}


void 
Compressible_PhaseField::SecondOrderSSPRK (Vector<MultiFab>& mf_new, 
                                           Vector<MultiFab>& mf_old, 
                                           int lev, 
                                           Real time, 
                                           Real dt_lev, 
                                           int ncycle)
{
    // swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }

#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }
    
    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif

    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);
    MultiFab flagfab    (grids[lev],dmap[lev],1,1);
    Real min_EnergyBound;



    // RK2 stage 1
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PreTimeStep(Uborder,ncycle,time,lev);

    compute_dUdt_FV(dUdt,Uborder,lev, time, Real(0.5)*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_new,ncycle,time+dt_lev,lev);


    /* RK2 stage 2 */
    // After fillpatch Uborder = U^n+dt*dUdt^n = U^*
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time+dt_lev,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    compute_dUdt_FV(dUdt,Uborder,lev, time+dt_lev, Real(0.5)*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    // U_new = 0.5*(Uborder+U_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(U_new, Real(0.5), Uborder, 0, Real(0.5), U_old, 0, 0, NSTATE, 0);
    // U_new += 0.5*dt*dUdt

    MultiFab::Saxpy(U_new, Real(0.5)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PostTimeStep(U_new,ncycle,time+dt_lev,lev);

}


void 
Compressible_PhaseField::ThirdOrderSSPRK (Vector<MultiFab>& mf_new, 
                                          Vector<MultiFab>& mf_old, 
                                          int lev, 
                                          Real time, 
                                          Real dt_lev, 
                                          int ncycle)
{

	// swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }



#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }
    
    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif




    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);
    // MultiFab flagfab    (grids[lev],dmap[lev],1,1);
    Real min_EnergyBound;



    // RK3 stage 1
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PreTimeStep(Uborder,ncycle,time,lev);


    compute_dUdt_FV(dUdt,Uborder,lev, time, (Real(1.0)/Real(6.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);
    
    MultiFab::LinComb(U_new, Real(1.0), Uborder, 0, dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }
    PostTimeStage(U_new,ncycle,time+dt_lev,lev);



    /* RK3 stage 2 */
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time+dt_lev,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);


    compute_dUdt_FV(dUdt,Uborder,lev, time+dt_lev, (Real(1.0)/Real(6.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    MultiFab::LinComb(U_new, Real(0.25), Uborder, 0, Real(0.75), U_old, 0, 0, NSTATE, 0);
    MultiFab::Saxpy(U_new, Real(0.25)*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }


    PostTimeStage(U_new,ncycle,time+Real(0.5)*dt_lev,lev);


    /* RK3 stage 3 */
    // U_new = U^(2) is at time+dt/2, not at t_new[lev] = time+dt, so fill it as a stage state
    stage = 2;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_new,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);


    compute_dUdt_FV(dUdt,Uborder,lev, time+Real(0.5)*dt_lev, (Real(2.0)/Real(3.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);
    
    MultiFab::LinComb(U_new, (Real(2.0)/Real(3.0)), Uborder, 0, (Real(1.0)/Real(3.0)), U_old, 0, 0, NSTATE, 0);
    MultiFab::Saxpy(U_new, (Real(2.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }


    PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PostTimeStep(U_new,ncycle,time+dt_lev,lev);
}


void 
Compressible_PhaseField::FourthOrderRK (Vector<MultiFab>& mf_new, 
                                        Vector<MultiFab>& mf_old, 
                                        int lev, 
                                        Real time, 
                                        Real dt_lev, 
                                        int ncycle)
{
	// swap the data, i.e. set dof_new = dof_old as we would for forward Euler
    std::swap(mf_old[lev], mf_new[lev]);

    int num_grow = NGROW;

   // Flux registers
   FluxRegister* fr_as_crse = nullptr;
   FluxRegister* fr_as_crse_nc = nullptr;
    if (do_reflux && lev < finest_level) {
        fr_as_crse = flux_reg[lev+1].get();
    }

    FluxRegister* fr_as_fine = nullptr;
    FluxRegister* fr_as_fine_nc = nullptr;
    if (do_reflux && lev > 0) {
        fr_as_fine = flux_reg[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }


#if (NONCONSERVATIVE == true)
   // Nonconservative flux registers
    if (do_reflux && lev < finest_level) {
        fr_as_crse_nc = flux_reg_nc[lev+1].get();
    }
    if (do_reflux && lev > 0) {
        fr_as_fine_nc = flux_reg_nc[lev].get();
    }

    if (fr_as_crse) {
        fr_as_crse_nc->setVal(Real(0.0));
    }
#endif



    // State
    MultiFab& U_new = mf_new[lev];
    MultiFab& U_old = mf_old[lev];
    MultiFab U_temp     (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab dUdt       (grids[lev],dmap[lev],NSTATE,0       );
    MultiFab Uborder    (grids[lev],dmap[lev],NSTATE,num_grow);
    MultiFab flagfab    (grids[lev],dmap[lev],1,1);
    Real min_EnergyBound;

    amrex::MultiFab::Copy(U_temp, U_new,0,0,NSTATE,0);


 
    //=====================================================================
    //RK4 stage 1 
    //=====================================================================
    int stage = 0;
    dUdt.setVal(Real(0.0));
    FillPatch(Uborder, lev, time,
              mf_new, t_new,
              mf_old, t_old,
              0, 0, NSTATE);

    PreTimeStep(Uborder,ncycle,time,lev);

    compute_dUdt_FV(dUdt,Uborder,lev, time, (Real(1.0)/Real(6.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    /* u_1 = u_n + 1/2 dt R(u_n) */
    MultiFab::LinComb(U_temp, Real(1.0), Uborder, 0, Real(0.5)*dt_lev, 			  dUdt, 0, 0, NSTATE, 0);
    /* u_* = u_n + 1/6 dt R(u_n) */
    MultiFab::LinComb(U_new,  Real(1.0), Uborder, 0, (Real(1.0)/Real(6.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_temp,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_temp,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }


    PostTimeStage(U_temp,ncycle,time+Real(0.5)*dt_lev,lev);


    //=====================================================================
    //RK4 stage 2 
    //=====================================================================
    // the stage state u_1 at time+dt/2 is in U_temp, while U_new holds the running sum u_*
    stage = 1;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    compute_dUdt_FV(dUdt,Uborder,lev, time+Real(0.5)*dt_lev, (Real(1.0)/Real(3.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    /* u_2 = u_n + 1/2 dt R(u_1) */
    MultiFab::LinComb(U_temp, Real(1.0), U_old, 0, Real(0.5)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    /* u_** = u_* + 1/3 dt R(u_1) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_temp,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_temp,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_temp,ncycle,time+Real(0.5)*dt_lev,lev);


 
    //=====================================================================
    //RK4 stage 3 
    //=====================================================================
    // the stage state u_2 at time+dt/2 is in U_temp
    stage = 2;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+Real(0.5)*dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    compute_dUdt_FV(dUdt,Uborder,lev, time+Real(0.5)*dt_lev, (Real(1.0)/Real(3.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    /* u_3 = u_n + dt R(u_2) */
    MultiFab::LinComb(U_temp, Real(1.0), U_old, 0, Real(1.0)*dt_lev, dUdt, 0, 0, NSTATE, 0);
    /* u_*** = u_** + 1/3 dt R(u_2) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(3.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_temp,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_temp,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_temp,ncycle,time+dt_lev,lev);



    //=====================================================================
    //RK4 stage 4
    //=====================================================================
    // the stage state u_3 at time+dt is in U_temp, not in U_new, so it is filled as a stage state
    stage = 3;
    dUdt.setVal(Real(0.0));
    FillPatchStage(Uborder, lev, time+dt_lev, U_temp,
                   mf_new, t_new,
                   mf_old, t_old,
                   0, 0, NSTATE);

    compute_dUdt_FV(dUdt,Uborder,lev, time+dt_lev, (Real(1.0)/Real(6.0))*dt_lev, ncycle, stage, fr_as_crse, fr_as_fine, fr_as_crse_nc, fr_as_fine_nc);

    /* u_n+1 = u_*** + 1/6 dt R(u_3) */
    MultiFab::Saxpy(U_new, (Real(1.0)/Real(6.0))*dt_lev, dUdt, 0, 0, NSTATE, 0);

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1){
        MinEnergyBound(min_EnergyBound,
                       U_new,
                       *h_parm);
        h_parm->FiniteVolume_Parm.ID_Return = 0;
        if (min_EnergyBound < 0){
            h_parm->FiniteVolume_Parm.ID_Return = 1;
            amrex::Real min_Pressure;
            MinPressure(min_Pressure,
                        U_new,
                        *h_parm);
            //amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<" < "<<d_parm->p_min<<"\n";
            amrex::Print()<<"\n"<<"min_EnergyBound from the updated state is "<<min_EnergyBound<<" < 0, and min_Pressure is "<<min_Pressure<<"\n";
            amrex::Print()<<"\n"<<"Restarting this time step"<<"\n";
            return;
        }
    }

    PostTimeStage(U_new,ncycle,time+dt_lev,lev);
    PostTimeStep(U_new,ncycle,time+dt_lev,lev);



}


void 
Compressible_PhaseField::AdvanceAtLevel (Vector<MultiFab>& mf_new, 
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

    if (h_parm->FiniteVolume_Parm.ID_Bound == 1 && h_parm->FiniteVolume_Parm.ID_Return == 1){
        return;
    }

}


void
Compressible_PhaseField::timeStepWithSubcycling (Vector<MultiFab>& mf_new,
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

                    // if there are newly created levels, set the time step. A level takes nsubsteps[k]
                    // steps per step of level k-1 (as in ComputeDt), so dt[k] = dt[k-1]/nsubsteps[k]
                    // brings it to the end time of level k-1 (a smaller dt would leave it behind)
                    for (int k = old_finest+1; k <= finest_level; ++k) {
                        dt[k] = dt[k-1] / nsubsteps[k];
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
    if (h_parm->FiniteVolume_Parm.ID_Bound == 1 && h_parm->FiniteVolume_Parm.ID_Return == 1){
        return;
    }

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
            if (h_parm->FiniteVolume_Parm.ID_Bound == 1 && h_parm->FiniteVolume_Parm.ID_Return == 1){
                return;
            }
        }

        if (do_reflux)
        {
            RefluxLev(lev);
        }

        AverageDownTo(mf_new,lev); // average lev+1 down to lev
    }
}


void 
Compressible_PhaseField::ExplicitTimeStep (amrex::Real cur_time)
{

    int lev = 0;
    int iteration = 1;
    int stage = 0;
    if (do_subcycle){
        timeStepWithSubcycling(dof_new, dof_old, lev, cur_time, iteration, stage);
        if (h_parm->FiniteVolume_Parm.ID_Bound == 1 && h_parm->FiniteVolume_Parm.ID_Return == 1){
            return;
        }
    }
    else{
        timeStepNoSubcycling(cur_time, iteration);
    }
}


void
Compressible_PhaseField::EstTimeStep (Real& dt_est, int lev, Real const& time)
{
    BL_PROFILE("Compressible_PhaseField::EstTimeStep()");

    dt_est = std::numeric_limits<Real>::max();

    const Real* dx  =  geom[lev].CellSize();

    if (time == 0.0 || h_parm->FiniteVolume_Parm.ID_Return == 1) {
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
}

void
Compressible_PhaseField::CalculateInitialcmax (int lev)
{

    const int num_grow = NGROW;

    Real time = 0.0;
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
    amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
    auto const prob_lo = Geom(lev).ProbLoArray();

    // State
    MultiFab& c_max_lev = c_max[lev];
    c_max_lev.setVal(0.0);

    // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], dof_new[lev].nComp(), num_grow);
    FillPatch(Sborder, lev, time,
              dof_new, t_new,
              dof_old, t_old,
              0, 0, Sborder.nComp());

    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(dof_new[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // ======== FLUX CALC AND UPDATE CMAX =========
            const Box& bx = mfi.tilebox();

            // Pull the data into an array
            Array4<Real const> statein   = Sborder.const_array(mfi);
            Array4<Real      > c_max_new = c_max_lev.array(mfi);
            amrex::ignore_unused(bx);

            for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                // the faces of this tile only: c_max(i,j,k) is written from the faces (i,j,k),
                // which are in the tile of cell (i,j,k), so no two OpenMP threads write one cell
                const Box bx_face = mfi.nodaltilebox(idim);
                amrex::ParallelFor(bx_face,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Array<Real,NSTATE> fhat;
                    Array<Real,NC_TERMS> fhatNC;
                    Real cmax;

                    FVM_Conservative2FluxHyperbolic_K(i, j, k,
                                                      fhat, fhatNC, cmax,
                                                      statein, *lparm,
                                                      dX, idim,
                                                      lparm->FiniteVolume_Parm, lparm->Physics_Parm);
                    
                    // amrex::Print() << cmax << "\n";
                    c_max_new(i,j,k,0) = std::max(cmax, c_max_new(i,j,k,0));

#if (DIFFUSION==true)
                    Array<Real,NSTATE> fhatD;
                    amrex::Real dmax;
{BL_PROFILE("compute_dUdt_FV::{ computing the idim-direction parabolic fluxes }");
#if (PHYSICS == SIXEQS)
                    FDM_Conservative2FluxDiffusion_K(i, j, k,
                                                     fhatD, fhatNC, dmax,
                                                     statein,
                                                     prob_lo,
                                        AMREX_D_DECL(dx,dy,dz),
                                                     time,
                                                     *lparm,
                                                     idim);
#else
                    FDM_Conservative2FluxDiffusion_K(i, j, k,
                                                     fhatD, dmax,
                                                     statein,
                                                     prob_lo,
                                        AMREX_D_DECL(dx,dy,dz),
                                                     time,
                                                     *lparm,
                                                     idim);
#endif
}
                    c_max_new(i,j,k,1) = std::max(dmax, c_max_new(i,j,k,1));
#endif

                });
            } // end idim
        } // end mfi
    } // end omp

    Print() << "Initial c_max on level " << lev << " is " << c_max[lev].norminf(0,0,false) << "\n";
#if (DIFFUSION==true)
    Print() << "Initial d_max on level " << lev << " is " << c_max[lev].norminf(1,0,false) << "\n";
#endif
}


