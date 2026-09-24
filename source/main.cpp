// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
// Portions derived from AMReX, BSD-3-Clause. See THIRD_PARTY_NOTICES.md.

#include <iostream>

#include <AMReX.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <COMPAS.H>
#include <Physics_DerivedClasses.H>
#include <Banner.H>


using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    PrintBanner();

    {

        // timer for profiling
        BL_PROFILE("main()");

        // wallclock time 
        const auto strt_total = amrex::second();

        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures

#if (PHYSICS == FIVEEQS)
        Compressible5Eq_PhaseField run;  

#elif (PHYSICS == FIVEEQS_NPHASE)
        Compressible5Eq_NPhase_PhaseField run;  

#elif (PHYSICS == SIXEQS)
        Compressible6Eq_PhaseField run;

#elif (PHYSICS == SIXEQS_IE_NPHASE)
        Compressible6Eq_IE_NPhase_PhaseField run;
#else
        COMPAS run;
#endif
        
        // initialize AMR data
        run.InitData();

        // advance solution to final time
        run.Evolve();

        const auto avg_timestep = run.avg_timestep;
        const auto evolve_total = run.evolve_total;


        // wallclock time
        auto end_total = amrex::second() - strt_total;

        if (run.Verbose()) {
            // print wallclock time
            ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "\n";
            if (run.timers == 1){
                amrex::Print() << "Average Timestep Time: " << avg_timestep << '\n';
                amrex::Print() << "Total Evolution Time: " << evolve_total << '\n';
            }
            amrex::Print() << "Total Time: " << end_total << '\n';
        }

    }

    amrex::Finalize();
}






