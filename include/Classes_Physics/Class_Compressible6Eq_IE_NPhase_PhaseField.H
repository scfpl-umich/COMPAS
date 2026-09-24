// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Compressible6Eq_IE_NPhase_PhaseField_H_
#define Compressible6Eq_IE_NPhase_PhaseField_H_

#include <COMPAS.H>
#include <Class_Compressible_PhaseField.H>

class Compressible6Eq_IE_NPhase_PhaseField
    : public Compressible_PhaseField
{
public:

    ~Compressible6Eq_IE_NPhase_PhaseField (); 

    void SetPhysicsBC () override;

    void TimeStepChecksDerived (int step) override;

    void PreTimeStep (MultiFab & mf, int step, amrex::Real time, int lev) override {};

    void PostTimeStage (MultiFab & mf, int step, amrex::Real time, int lev) override;

    void PostTimeStep (MultiFab & mf, int step, amrex::Real time, int lev) override {};

    void RefluxLev (int lev) override;

    void PhaseField_PostTimeStage (MultiFab & mf, int step, amrex::Real time, int lev) override;

    void PhaseField_RefluxLev (int lev) override;

    void SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab, amrex::MultiFab const & dof) override;

    void RefluxSourceTerms (amrex::MultiFab & new_dof,
                            amrex::MultiFab & source_fab,
                            int lev,
                            amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr) override;

    void PhaseField_SaveSourceTermsBeforeReflux (amrex::MultiFab & source_fab, 
                                                 amrex::MultiFab const & dof) override;

    void PhaseField_RefluxSourceTerms (amrex::MultiFab & new_dof,
                            amrex::MultiFab const & source_fab,
                            int lev,
                            amrex::Vector<std::unique_ptr<amrex::FluxRegister>> & fr) override;

    // =============================================================================================
    // Original system 
    // =============================================================================================

    void
    Correction_Energy (MultiFab& mf,
			           Parm const* parm);

    void
    Relaxation_Pressure (MultiFab& mf,
				         Parm const* parm);

};


// #endif
#endif