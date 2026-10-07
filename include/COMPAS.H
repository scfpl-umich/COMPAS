// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
// Portions derived from AMReX, BSD-3-Clause. See THIRD_PARTY_NOTICES.md.


#ifndef COMPAS_H_
#define COMPAS_H_

#include <string>
#include <limits>
#include <memory>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif


#include<AMReX_Gpu.H>
#include <AMReX_GpuMemory.H>
#include<AMReX_GpuComplex.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>
#include <AMReX_TimeIntegrator.H>
#include <AMReX_Vector.H>
#include <AMReX_Print.H>
#include <AMReX_FillPatchUtil.H>

/* Include custom packages */
#include <Macros.H>
#include <Tools.H>
#include <Parm.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <IncludePhysics.H>
#include <IncludeOperators.H>
#include <Tagging.H>
#include <InitData.H>
#include <LinearAlgebra.H>

#include <PreConstructors.H>


struct base_physics {};

struct OutputWindow
{
    amrex::Real t_lo;
    amrex::Real t_hi;
    amrex::Real dt;
};


template <typename F>
AMREX_FORCE_INLINE
void device_modify(F&& f)
{
    #if AMREX_USE_GPU == true
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
    {
        f();
    });
    #endif
}


class COMPAS
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    COMPAS ();
    virtual ~COMPAS();

    virtual void Physics_Derived_Constructor ();

    virtual void ReadParameters_Derived (){};

    virtual void MakeNewLevelFromCoarse_Derived (int lev, 
                                                 Real time, 
                                                 const BoxArray& ba,
                                                 const DistributionMapping& dm){};

    virtual void MakeNewLevelFromScratch_Derived (int lev, 
                                                  Real time, 
                                                  const BoxArray& ba,
                                                  const DistributionMapping& dm){};
    virtual void RemakeLevel_Derived (int lev, 
                                      Real time, 
                                      const BoxArray& ba,
                                      const DistributionMapping& dm){};

    virtual void RemakeAllLevelsFromData_Derived (int lev,
                                                  const BoxArray& ba,
                                                  const DistributionMapping& dm){};

    virtual void ClearLevel_Derived (int lev){};

    virtual void SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab,
                                              amrex::MultiFab const & dof) {};

    virtual void RefluxSourceTerms (amrex::MultiFab & new_dof,
                                    amrex::MultiFab & source_fab,
                                    int lev,
                                    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr) {};


    // advance solution to final time
    virtual void Evolve ();

    // initializes multilevel data
    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Remake all levels using provided input data
    void RemakeAllLevelsFromData (Vector<MultiFab> const& dof_Input,
                                  Vector<Real>     const&   t_Input,
                                  Vector<Real>     const&  dt_Input,
                                  Vector<int>      const& istep_Input,
                                  int              finest_level_Input);

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

     // Advance at a single level for a single time step, update flux registers
    virtual void AdvanceAtLevel (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle);

    virtual
    void ForwardEuler (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle) {};

    virtual
    void SecondOrderSSPRK (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle) {};

    virtual
    void ThirdOrderSSPRK (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle) {};

    virtual
    void FourthOrderRK (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle) {};

    // time is the time of the stage state mf_old, which the user source term receives
    virtual void
    compute_dUdt_FV (MultiFab& mf_new, 
                     MultiFab& mf_old, 
                     int lev, 
                     Real time,
                     Real dt_lev, 
                     int ncycle, 
                     int stage,
                     FluxRegister* fr_as_crse, 
                     FluxRegister* fr_as_fine,
                     FluxRegister* fr_as_crse_nc, 
                     FluxRegister* fr_as_fine_nc) {};

    // Advance at all levels for a single time step (the step without subcycling)
    virtual void AdvanceAllLevels (amrex::Real time, amrex::Real dt_lev, int iteration);

    // compute dt from CFL considerations
    virtual
    void EstTimeStep (amrex::Real& dt_est, int lev, amrex::Real const& time);


    // Initialize pointer and device variables
    void variableSetUp ();

    // Delete pointer and device variables
    void variableCleanUp ();

    void CalculateOutputMF ();

    void ResetOutputMF ();

    void OutputError ();

    void EvaluateError (amrex::Vector<amrex::Real> & L2, 
                        amrex::Vector<amrex::Real> & Linf);

    virtual void SetPhysicsBC (){};

    virtual void ExplicitTimeStep (amrex::Real cur_time);

    void TimeStepChecks (int step);

    virtual void TimeStepChecksDerived (int step){};

    virtual void PreTimeStep (MultiFab & mf, int step, amrex::Real time, int lev){};

    virtual void PostTimeStage (MultiFab & mf, int step, amrex::Real time, int lev){};

    virtual void PostTimeStep (MultiFab & mf, int step, amrex::Real time, int lev){};

    virtual void RefluxLev (int lev);

    void SaveLine (amrex::Vector<amrex::MultiFab> & mf, 
                   int dir, 
                   amrex::Array<amrex::Real,2> dir_coord, 
                   amrex::Vector<std::string> vnames,
                   std::string fname);

    void SavePlane (amrex::Vector<amrex::MultiFab> & mf, 
                    int dir, 
                    amrex::Real dir_coord, 
                    amrex::Vector<std::string> vnames,
                    std::string fname);

    amrex::Real avg_timestep = 0.0;
    amrex::Real evolve_total = 0.0;
    amrex::Real evolve_step = 0.0;

     // Choice of timestep output style - not connected with other verbosity parameters
    int timers = 1;

    // Vector for rolling average
    amrex::Vector<amrex::Real> ra_estimates;

    // Interval for rolling average
    int run_time_ra_int = 10;


#if !defined(AMREX_USE_CUDA)
// protected:
#endif

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    virtual void ReadParameters();

    void Coarsen (int ref_ratio);

    void Coarsen (int ref_ratio, int lev);

    void Refine (int ref_ratio);

    void Refine (int ref_ratio, int lev);

    //For cell-center data: average down data to crse_lev
    void AverageDownTo (Vector<MultiFab>& mf, int crse_lev);

    //For cell-center data: set covered coarse cells to be the average of overlying fine cells for all levels
    void AverageDown (Vector<MultiFab>& mf);

    //For face-center data: average down data to crse_lev
    void AverageDownFacesTo (Vector< Array<MultiFab,AMREX_SPACEDIM> >& mf, int crse_lev);

    //For face-center data: set covered coarse cells to be the average of overlying fine cells for all levels
    void AverageDownFaces (Vector< Array<MultiFab,AMREX_SPACEDIM> >& mf);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (MultiFab& mf,
                    int lev,
                    Real const& time,
                    Vector<MultiFab>& mf_new,
                    Vector<Real>&      t_new,
                    Vector<MultiFab>& mf_old,
                    Vector<Real>&      t_old,
                    int scomp,
                    int dcomp,
                    int ncomp);

    void FillPatch (MultiFab& mf,
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
                    Interpolater* mapper);

    // fill mf for a Runge-Kutta stage that starts at stage_time from the stage state mf_stage
    // (valid region and same-level ghost cells), with the physical boundary conditions at
    // stage_time and the coarser level interpolated in time to stage_time
    void FillPatchStage (MultiFab& mf,
                         int lev,
                         Real const& stage_time,
                         MultiFab& mf_stage,
                         Vector<MultiFab>& mf_new,
                         Vector<Real>&      t_new,
                         Vector<MultiFab>& mf_old,
                         Vector<Real>&      t_old,
                         int scomp,
                         int dcomp,
                         int ncomp);

    // fill mf of level lev for a Runge-Kutta stage that every level takes at the same time (no
    // subcycling): the valid region and same-level ghost cells from U_stage[lev], the physical
    // boundary conditions at stage_time, and the coarse-fine ghost cells interpolated in space from
    // the coarser level's stage state U_stage[lev-1] (no interpolation in time)
    void FillPatchStageLevels (MultiFab& mf,
                               int lev,
                               Real const& stage_time,
                               Vector<MultiFab>& U_stage,
                               int scomp,
                               int dcomp,
                               int ncomp);


    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (MultiFab& mf,
                          int lev,
                          Real const& time,
                          Vector<MultiFab>& mf_new,
                          Vector<Real>&      t_new,
                          Vector<MultiFab>& mf_old,
                          Vector<Real>&      t_old,
                          int scomp,
                          int dcomp,
                          int ncomp);

    void FillCoarsePatch (MultiFab& mf,
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
                          Interpolater* mapper);

    // utility to copy in data from dof_old and/or dof_new into another multifab
    void GetData (Vector<MultiFab*>& data,
                  Vector<Real>& datatime,
                  int lev,
                  Real const& time,
                  Vector<MultiFab>& mf_new,
                  Vector<Real>&      t_new,
                  Vector<MultiFab>& mf_old,
                  Vector<Real>&      t_old);

    // Advance a level by dt - includes a recursive call for finer levels
    virtual
    void timeStepWithSubcycling (Vector<MultiFab>& mf_new,
                                 Vector<MultiFab>& mf_old,
                                 int lev, 
                                 amrex::Real time, 
                                 int iteration,
                                 int stage);

    // Advance all levels by the same dt
    void timeStepNoSubcycling (amrex::Real time, int iteration);

    // Set the time step to a user-specified time step (for convergence testing)
    void FixedDt ();

    // a wrapper for EstTimeStep
    void ComputeDt ();

    // Get the initial wave speeds
    virtual void CalculateInitialcmax (int lev) {};


    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // void get_data_at_points (Vector<MultiFab>& mf_vec);

    // void get_data_at_lines (Vector<MultiFab>& mf_vec);

    // void get_data_at_planes (Vector<MultiFab>& mf_vec);

    bool OutputTimeInterval (amrex::Real cur_time);

    void SumVar ();

    void TimestepOutput (int step, amrex::Real cur_time);

    void ResolvePaths();

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxFinest();


    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> dof_new;
    amrex::Vector<amrex::MultiFab> dof_old;

    int ID_IB = 0;
    amrex::Vector<amrex::MultiFab> ib_levelset;


    amrex::Vector<amrex::MultiFab>* mf_new;
    amrex::Vector<amrex::MultiFab>* mf_old;


    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> bcs; 

    amrex::Vector<amrex::BCRec> flag_bcs; 

    amrex::Vector<int> lo_bc;
    amrex::Vector<int> hi_bc;

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> >    flux_reg;
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> >    flux_reg_nc;

    // Velocity on all faces at all levels
    // amrex::Vector< amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > facevel;

    // Maximum wave speed on each level
    amrex::Vector<amrex::MultiFab> c_max;

    // Multifab for storing element-local flags 
    amrex::Vector<amrex::MultiFab> flags;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // Choose the coordinate system, default is Cartesian
    int ID_COORDSYS = 0;

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.4;
    amrex::Real cfl_fast = false;
    int cfl_switch = false;

    // diffusive vnn number - dt = vnn*dx*dx/dmax
    amrex::Real vnn = 0.25;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // do we subcycle in time?
    int do_subcycle = 1;

    // plotfile prefix and frequency
    std::string plot_file = "";
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

    // UserOutputFunction is called every user_output_int coarse steps,
    // or with each plotfile when user_output_int <= 0
    int user_output_int = 0;

    bool ID_QuadratureIC = true; // Use quadrature for the initial condition?


    // Vector of TagVar structs that each contain information about refinement for a given variable
    amrex::Vector<TagVar> tag_vars;

    // Vector of variables to be output
    amrex::Vector<std::string> output_vars;

    // Number of output variables
    int num_output_vars = NSTATE;

    // Boolean for outputing the state variables
    bool output_state_vars = true;

    // Array of indices for output variables
    amrex::Array<int, 2*NSTATE> out_vars;
    amrex::Vector<std::string> out_vars_string = {};

    //total number of components
    int nC;

    // Index of component to be summed (conservation check)
    int sum_var = 0;

    // Checks the cfl estimate in the residual calculation
    int check_cfl = 1;

    // Specific output times (within a dt interval)
    amrex::Vector<amrex::Real> t_write;

    // Interval for outputs (within a dt interval)
    bool t_write_interval_true = false;
    amrex::Real t_write_interval = 0.0;

    amrex::Vector<OutputWindow> output_windows;

    // User-specified dt (for convergence testing)
    amrex::Real dt_fixed = 0.0;

    // Name of output file for convergence error data 
    std::string conv_output_file = "err";

    // Name of the test problem for convergence
    std::string conv_prob_name;

    int output_subdirectory = 0;
    std::string output_dir = "";
    std::string case_name = "";
    std::string movie_filename = "";

    // Start and stop times for the simulation
    amrex::Real Evolve_start_time;
    amrex::Real Evolve_stop_time;

    // Whether or not to sum a variable (costly, so pay attention for large sims)
    int timestep_sum = 1;

    // The member to store the sum 
    amrex::Real sum_state;

    std::string ID_TimeIntegrator = "TVD-RK3"; // Default to SSP-RK3

    std::string ID_Interpolater = "lincc_interp";

    // Whether or not to limit the change in dt between time steps
    int timestep_change_limiter = 1;

    // Number of cells in each direction
    amrex::Vector<int> n_cells; 

    // Number of points to save data at (closest cell-center)
    int num_save_points = 0;

    // Interval for saving points 
    int save_points_int = 0;

    // Points to save
    amrex::Vector<amrex::Real> points_to_save;

    // File prefix for each point being saved
    std::string points_filename = "";

    // Number of straight lines to save data at (closest cell-center)
    int num_save_straight_lines = 0;

    // x-y point pairs for a z-line
    amrex::Vector<amrex::Real> z_lines_to_save;

    std::string lines_filename = "";

    int save_lines_int = 0;

    base_physics physics;

public:
    static Parm* h_parm;
    static Parm* d_parm;

    int NGROW;

};



#endif
