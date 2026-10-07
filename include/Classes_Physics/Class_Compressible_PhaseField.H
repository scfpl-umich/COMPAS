// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Compressible_PhaseField_H_
#define Compressible_PhaseField_H_

// #include<Macros.H>
// #if (PHYSICS==FIVEEQS)

#include <COMPAS.H>

class Compressible_PhaseField
    : public COMPAS
{
public:

    
    ~Compressible_PhaseField (); 

    void Evolve () override;

    void Physics_Derived_Constructor () override;

    void ReadParameters_Derived () override;

    void MakeNewLevelFromCoarse_Derived (int lev, 
                                         Real time, 
                                         const BoxArray& ba,
                                         const DistributionMapping& dm) override;

    void MakeNewLevelFromScratch_Derived (int lev, 
                                          Real time, 
                                          const BoxArray& ba,
                                          const DistributionMapping& dm) override;
    void RemakeLevel_Derived (int lev, 
                              Real time, 
                              const BoxArray& ba,
                              const DistributionMapping& dm) override;

    void RemakeAllLevelsFromData_Derived (int lev,
                                          const BoxArray& ba,
                                          const DistributionMapping& dm) override;

    void ClearLevel_Derived (int lev) override;

    void SetPhysicsBC () override;

    void CalculateInitialcmax (int lev) override;

    // =============================================================================================
    // Original system 
    // =============================================================================================
    void ExplicitTimeStep (amrex::Real cur_time) override;

    void timeStepWithSubcycling (Vector<MultiFab>& mf_new,
                                 Vector<MultiFab>& mf_old,
                                 int lev, 
                                 amrex::Real time, 
                                 int iteration,
                                 int stage) override;

    void AdvanceAtLevel (Vector<MultiFab>& mf_new, 
                         Vector<MultiFab>& mf_old, 
                         int lev, 
                         Real time, 
                         Real dt_lev, 
                         int ncycle) override;

    void ForwardEuler (Vector<MultiFab>& mf_new, 
                       Vector<MultiFab>& mf_old, 
                       int lev, 
                       Real time, 
                       Real dt_lev, 
                       int ncycle) override;

    void SecondOrderSSPRK (Vector<MultiFab>& mf_new, 
                           Vector<MultiFab>& mf_old, 
                           int lev, 
                           Real time, 
                           Real dt_lev, 
                           int ncycle) override;

    void ThirdOrderSSPRK (Vector<MultiFab>& mf_new, 
                          Vector<MultiFab>& mf_old, 
                          int lev, 
                          Real time, 
                          Real dt_lev, 
                          int ncycle) override;

    void FourthOrderRK (Vector<MultiFab>& mf_new, 
                        Vector<MultiFab>& mf_old, 
                        int lev, 
                        Real time, 
                        Real dt_lev, 
                        int ncycle) override;

    void
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
                     FluxRegister* fr_as_fine_nc) override;

    // Parts of compute_dUdt_FV, also used by the step without subcycling: the face MultiFabs of a
    // level, the face fluxes and NC face quantities (and c_max) from a state with ghost cells, the
    // right-hand side from the face MultiFabs, the CFL check, and the area scaling for the registers
    void DefineFaceFluxes (Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                           Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                           int lev);

    void ComputeFaceFluxes (MultiFab const& U,
                            Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                            Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                            int lev,
                            Real time,
                            int stage);

    void FluxDivergence (MultiFab& dUdt_mf,
                         MultiFab const& U,
                         Array<MultiFab,AMREX_SPACEDIM> const& fluxes,
                         Array<MultiFab,AMREX_SPACEDIM> const& fluxes_nc,
                         int lev,
                         Real time,
                         Real dt_lev);

    void CheckCFL (int lev, Real dt_lev);

    void ScaleFaceFluxesByArea (Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                int lev);

    // =============================================================================================
    // Step without subcycling, all levels stage by stage (Class_Compressible_PhaseField_StageCoupling.cpp)
    // =============================================================================================
    void AdvanceAllLevels (Real time, Real dt_lev, int ncycle) override;

    // One step of all levels: Grad_Q == nullptr for the hyperbolic step, else the Phase-Field step
    void AdvanceAllLevelsStageCoupled (Vector<MultiFab>& mf_new,
                                       Vector<MultiFab>& mf_old,
                                       Vector<Array<MultiFab,AMREX_SPACEDIM>> const* Grad_Q,
                                       Real time,
                                       Real dt_lev,
                                       int ncycle);

    // Right-hand side of every level for one stage, with the finer faces averaged onto the coarser
    void StageRHSAllLevels (Vector<MultiFab>& dUdt,
                            Vector<MultiFab>& Uborder,
                            Vector<MultiFab>& U_stage,
                            Vector<Array<MultiFab,AMREX_SPACEDIM>> const* Grad_Q,
                            Real stage_time,
                            Real dt_stage,
                            int ncycle,
                            int stage);

    bool StageBoundViolated (MultiFab const& U);

    void EstTimeStep (Real& dt_est, int lev, Real const& time) override;

    virtual
    void MinEnergyBound (amrex::Real & min_EnergyBound,
                         amrex::MultiFab const& Conservative,
                         Parm const& EOS_Parameter);

    virtual
    void MinEnergyBound (amrex::Real & min_EnergyBound,
                         amrex::Vector<amrex::MultiFab> const& Conservative,
                         Parm const& EOS_Parameter);

    virtual
    void MinPressure (amrex::Real & min_Pressure,
                      amrex::MultiFab const& Conservative,
                      Parm const& EOS_Parameter);

    virtual
    void MinPressure (amrex::Real & min_Pressure,
                      amrex::Vector<amrex::MultiFab> const& Conservative,
                      Parm const& EOS_Parameter);

    virtual
    void MaxVelocity (amrex::Real & max_V,
                      amrex::MultiFab const& Conservative);

    virtual
    void MaxVelocity (amrex::Real & max_V,
                      amrex::Vector<amrex::MultiFab> const& Conservative);
    // =============================================================================================
    // Phase-Field
    // =============================================================================================

    // Compute Div_J from the conservative Allen-Cahn model
    virtual
    void PhaseField_Mechanism (Vector<MultiFab>& Div_J,
                               Vector<MultiFab> const& PFBorder,
                               Vector<MultiFab> const& PsiBorder,
                               PhaseField_Parameter const& PhaseField_Parm);

    virtual
    void PhaseField_Mechanism (Vector<MultiFab>& Div_J,
                               Vector<Array<MultiFab,AMREX_SPACEDIM>>& J_Flux,
                               Vector<MultiFab> const& PFBorder,
                               Vector<MultiFab> const& PsiBorder,
                               PhaseField_Parameter const& PhaseField_Parm);

    // Compute J_Flux from the Cahn-Hilliard model
    virtual
    void PhaseField_CH (Vector<Array<MultiFab,AMREX_SPACEDIM>>& J_Flux,
                        Vector<MultiFab> const& PFBorder,
                        Vector<MultiFab> const& PsiBorder,
                        Real const& Mobility,
                        Real const& Lambda,
                        Real const& Eta,
                        int ID_Mechanism,
                        int ID_Potential);

    // Compute J_Flux from the conservative Diffuse-Interface (Phase-Field) model
    virtual
    void PhaseField_CDI (Vector<Array<MultiFab,AMREX_SPACEDIM>>& J_Flux,
                         Vector<MultiFab> const& PFBorder,
                         Vector<MultiFab> const& PsiBorder,
                         Real const& Mobility,
                         Real const& Eta,
                         int ID_Mechanism);

    // Compute J_Flux from the accurate conservative Diffuse-Interface model
    virtual
    void PhaseField_ACDI (Vector<Array<MultiFab,AMREX_SPACEDIM>>& J_Flux,
                          Vector<MultiFab> const& PFBorder,
                          Vector<MultiFab> const& PsiBorder,
                          Real const& Mobility,
                          Real const& Eta,
                          int ID_Mechanism);

    // Compute diffused Heaviside function based on the Phase-Field function
    virtual
    void PhaseField_Diffusion (Vector<MultiFab>& HeavisideBorder,
                               Vector<MultiFab> const& PF,
                               Real const& Eta);

    // Compute signed distance function based on the Phase-Field function
    virtual
    void PhaseField_Reinitialization (Vector<MultiFab>& LSBorder,
                                      Vector<MultiFab> const& PF,
                                      Real const& Eta);

    // Compute Normal vector at cell centers from function f
    virtual
    void PhaseField_Normal (Vector<MultiFab>& NormalBorder,
                            Vector<MultiFab> const& fBorder,
                            int ID_Phase);

    // Compute Grad_Q from Div_J
    virtual
    void PhaseField_GradQPF (Vector<Array<MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                             Vector<MultiFab> const& Div_J,
                             Vector<MultiFab> const& PFBorder,
                             amrex::Real const& Eta,
                             amrex::Real const& Eta_Factor,
                             int ID_Weight,
                             LinearSystem_Parameter const& LinearSystem_Parm);

    // Compute Grad_QPF from J_Flux
    virtual
    void PhaseField_GradQPF (Vector<Array<MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                             Vector<Array<MultiFab,AMREX_SPACEDIM>> const& J_Flux,
                             Vector<MultiFab> const& PFBorder,
                             Real const& Eta,
                             Real const& Eta_Factor);

    // Compute Grad_QPF from J_Flux with degenerate mobility
    virtual
    void PhaseField_GradQPF (Vector<Array<MultiFab,AMREX_SPACEDIM>>& Grad_QPF,
                             Vector<Array<MultiFab,AMREX_SPACEDIM>> const& J_Flux);

    //Advance state vector at a level with the Phase-Field flux vector
    virtual
    void PhaseField_compute_dUdt_FV (MultiFab& mf_new, 
                                     MultiFab& mf_old, 
                                     Array<MultiFab,AMREX_SPACEDIM> const& Grad_Q, 
                                     int lev, 
                                     Real dt_lev, 
                                     int ncycle, 
                                     int stage,
                                     FluxRegister* fr_as_crse, 
                                     FluxRegister* fr_as_fine, 
                                     FluxRegister* fr_as_crse_nc, 
                                     FluxRegister* fr_as_fine_nc);

    // Parts of PhaseField_compute_dUdt_FV: the Phase-Field face fluxes of a level, and the
    // right-hand side from them
    void PhaseField_ComputeFaceFluxes (MultiFab const& U,
                                       Array<MultiFab,AMREX_SPACEDIM> const& Grad_Q,
                                       Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                                       Array<MultiFab,AMREX_SPACEDIM>& fluxes_nc,
                                       int lev,
                                       int stage);

    void PhaseField_FluxDivergence (MultiFab& dUdt_mf,
                                    MultiFab const& U,
                                    Array<MultiFab,AMREX_SPACEDIM> const& fluxes,
                                    Array<MultiFab,AMREX_SPACEDIM> const& fluxes_nc,
                                    int lev);

    //Euler advancement
    void PhaseField_ForwardEuler (Vector<MultiFab>& mf_new, 
                                   Vector<MultiFab>& mf_old, 
                                   Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                   int lev, 
                                   Real time, 
                                   Real dt_lev, 
                                   int ncycle);

    //RK2 advancement
    void PhaseField_SecondOrderSSPRK (Vector<MultiFab>& mf_new, 
                                      Vector<MultiFab>& mf_old, 
                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                      int lev, 
                                      Real time, 
                                      Real dt_lev, 
                                      int ncycle);

    //RK3 advancement
    void PhaseField_ThirdOrderSSPRK (Vector<MultiFab>& mf_new, 
                                     Vector<MultiFab>& mf_old, 
                                     Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                     int lev, 
                                     Real time, 
                                     Real dt_lev, 
                                     int ncycle);

    //RK4 advancement
    void PhaseField_FourthOrderRK (Vector<MultiFab>& mf_new, 
                                   Vector<MultiFab>& mf_old, 
                                   Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                   int lev, 
                                   Real time, 
                                   Real dt_lev, 
                                   int ncycle);

    virtual void PhaseField_PreTimeStep (MultiFab & mf, int step, amrex::Real time, int lev){};

    virtual void PhaseField_PostTimeStage (MultiFab & mf, int step, amrex::Real time, int lev){};

    virtual void PhaseField_PostTimeStep (MultiFab & mf, int step, amrex::Real time, int lev){};

    //Advance state vector at a level with different RK schemes
    void PhaseField_AdvanceAtLevel (Vector<MultiFab>& mf_new, 
                                    Vector<MultiFab>& mf_old, 
                                    Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                                    int lev, 
                                    Real time, 
                                    Real dt_lev, 
                                    int ncycle);

    //Advance state vector at all levels, stage by stage (the step without subcycling)
    void PhaseField_AdvanceAllLevels (Vector<MultiFab>& mf_new,
                                      Vector<MultiFab>& mf_old,
                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                      Real time,
                                      Real dt_lev,
                                      int ncycle);

    // Advance a level by PF_dt (includes a recursive call for finer levels)
    void PhaseField_TimeStepWithSubcycling (Vector<MultiFab>& mf_new,
                                            Vector<MultiFab>& mf_old,
                                            Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                            int lev, 
                                            amrex::Real time, 
                                            int iteration,
                                            int stage);

    // Advance all the levels with the same PF_dt
    void PhaseField_TimeStepNoSubcycling (Vector<MultiFab>& mf_new,
                                          Vector<MultiFab>& mf_old,
                                          Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q,
                                          Real time,
                                          int iteration);

    virtual void PhaseField_RefluxLev (int lev);

    //Explicit advancement in a time step
    void PhaseField_ExplicitTimeStep (amrex::Real cur_time,
                                      Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q);

    //Compute Phase-Field time step
    virtual
    void PhaseField_ComputeDt (Vector<Array<MultiFab,AMREX_SPACEDIM>> const& Grad_Q, 
                               Real PF_time);

    // advance solution to final time
    virtual
    void PhaseField_Evolve ();


    // advance Phase-Field mechanism to one hyperbolic time step
    virtual
    void PhaseFieldMechanism_Evolve (Vector<MultiFab>& Div_Jhat,
                                     Vector<Array<MultiFab,AMREX_SPACEDIM>>& Jhat_Flux,
                                     Vector<MultiFab> const& PF,
                                     Vector<BCRec> const& PF_bcs,
                                     PhaseField_Parameter const& PhaseField_Parm);


    virtual
    void PhaseField_SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab, 
                                                 amrex::MultiFab const & dof) {};

    virtual 
    void PhaseField_RefluxSourceTerms (amrex::MultiFab & new_dof,
                            amrex::MultiFab const & source_fab,
                            int lev,
                            amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr) {};

    // =============================================================================================
    //                                  Bound-preserving for the Hyperbolic step
    // =============================================================================================
/*
    int ID_Bound = 0; //0: Bound-preserving is off; 1: Bound-preserving is on
    int ID_Return = 0;
    int ID_Lambda = 1; //0: General; 1: Special (like all gases or single-phase flows)
    amrex::Real Quadrature_Weight = 1.0/6.0; //cmax*dt/dx/Quadrature_Weight=CFL
*/
    
    // =============================================================================================
    //                                  Phase-Field parameters
    // =============================================================================================
/*
    Real PF_Mob = d_parm->PhaseField_Parm.Mob;
*/

    amrex::Vector<amrex::Real> PF_dt;

    // Maximum wave speed on each level for the phase field 
    amrex::Vector<amrex::MultiFab> PF_c_max;


    // =============================================================================================
    //                                  Input deck parameters
    // =============================================================================================

};





// #endif
#endif