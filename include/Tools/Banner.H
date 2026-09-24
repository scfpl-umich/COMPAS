// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef BANNER_H_
#define BANNER_H_

// =======================================
// The banner printed at the start of a run
//
// PrintBanner() writes the name, the build that was compiled (model, dimension, parallel
// backends), the AMReX version and the rank count, and the license notice. It prints on the
// I/O rank only, and run.banner = 0 in the inputs file turns it off.
// =======================================

#include <string>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Version.H>
#include <Macros.H>

inline void PrintBanner()
{
    amrex::ParmParse pp("run");
    int banner = 1;
    pp.query("banner", banner);
    if (banner == 0) { return; }

    // ---------------------------------------
    // The model and the backends this executable was compiled with
    // ---------------------------------------
#if   (PHYSICS == FIVEEQS)
    std::string model = "FIVEEQS, five-equation, two phases";
#elif (PHYSICS == FIVEEQS_NPHASE)
    std::string model = "FIVEEQS_NPHASE, five-equation, " + std::to_string(NPHASE) + " phases";
#elif (PHYSICS == SIXEQS)
    std::string model = "SIXEQS, six-equation, two phases";
#elif (PHYSICS == SIXEQS_IE_NPHASE)
    std::string model = "SIXEQS_IE_NPHASE, six-equation, " + std::to_string(NPHASE) + " phases";
#else
    std::string model = "unknown";
#endif

    std::string parallel = "MPI x" + std::to_string(amrex::ParallelDescriptor::NProcs());
#ifdef AMREX_USE_OMP
    parallel += " + OpenMP";
#endif
#if   defined(AMREX_USE_CUDA)
    parallel += " + CUDA";
#elif defined(AMREX_USE_HIP)
    parallel += " + HIP";
#elif defined(AMREX_USE_SYCL)
    parallel += " + SYCL";
#endif

    // ---------------------------------------
    // The banner
    // ---------------------------------------
    amrex::Print()
        << "\n"
        << R"(     __________  __  _______  ___   _____)"         << "\n"
        << R"(    / ____/ __ \/  |/  / __ \/   | / ___/)"         << "\n"
        << R"(   / /   / / / / /|_/ / /_/ / /| | \__ \)"          << "\n"
        << R"(  / /___/ /_/ / /  / / ____/ ___ |___/ /)"          << "\n"
        << R"(  \____/\____/_/  /_/_/   /_/  |_/____/)"           << "\n"
        << "\n"
        << "  COmpressible Multiphase multiPhysics Adaptive Solver\n"
        << "  Scientific Computing and Flow Physics Laboratory, University of Michigan\n"
        << "\n"
        << "  Model     " << model << "\n"
        << "  Build     " << AMREX_SPACEDIM << "D, " << parallel << ", AMReX " << amrex::Version() << "\n"
        << "\n"
        << "  Copyright (C) 2026 The Regents of the University of Michigan.\n"
        << "  Free software under the GNU General Public License, version 3 or later, with\n"
        << "  ABSOLUTELY NO WARRANTY. See LICENSE for the terms, and README.md for the\n"
        << "  papers to cite.\n"
        << "\n";
}

#endif  // BANNER_H_
