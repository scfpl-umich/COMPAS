// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef LOG_UTILS_H_
#define LOG_UTILS_H_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <vector>
#include <utility>
#include <sys/stat.h>   // mkdir

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Version.H>

// -------------------------------------------------------------
// Simple path utilities (no <filesystem>)
// -------------------------------------------------------------
inline std::string PathDirname(const std::string& path)
{
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) return std::string{};
    if (pos == 0) return std::string{"/"};
    return path.substr(0, pos);
}

inline std::string PathBasename(const std::string& path)
{
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) return path;
    return path.substr(pos + 1);
}

inline std::string PathJoin(const std::string& a, const std::string& b)
{
    if (a.empty() || a == ".") return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

// Recursively create directories (mkdir -p behavior)
inline void CreateDirectories(const std::string& path)
{
    if (path.empty() || path == ".")
        return;

    std::string current;
    std::stringstream ss(path);
    std::string item;

    if (path[0] == '/') current = "/";

    while (std::getline(ss, item, '/')) {
        if (item.empty()) continue;
        current = PathJoin(current, item);
        mkdir(current.c_str(), 0755); // ignore errors
    }
}

// -------------------------------------------------------------
// Get timestamp
// -------------------------------------------------------------
inline std::string CurrentTimeString()
{
    std::time_t t = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return buf;
}

// -------------------------------------------------------------
// Safe loader: returns empty string if file doesn't exist.
// -------------------------------------------------------------
inline
std::string SafeLoadFileToString(const std::string& fname)
{
    std::ifstream f(fname, std::ios::in | std::ios::binary);
    if (!f.is_open())
        return {};

    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// -------------------------------------------------------------
// Runtime parameter dump (AMReX ParmParse)
// -------------------------------------------------------------
inline
std::string DumpRuntimeParameters()
{
    std::ostringstream ss;
    amrex::ParmParse::dumpTable(ss);
    return ss.str();
}

// -------------------------------------------------------------
// Main logging function (rank 0 only).
// Sections: { "Title", "path/to/file" }
// Missing files are silently skipped.
// -------------------------------------------------------------
inline
void DumpConfigurationLog(
    const std::string& plot_file,
    const std::string& case_name,
    const std::vector<std::pair<std::string,std::string>>& sections)
{
    using namespace amrex;

    if (ParallelDescriptor::MyProc() != 0)
        return;

    // ---- Determine log directory from plot_file ----
    std::string plot_dir = PathDirname(plot_file);
    if (plot_dir.empty()) plot_dir = ".";

    CreateDirectories(plot_dir);

    std::string log_path =
        PathJoin(plot_dir, "config_log_" + case_name + ".txt");

    std::ofstream log(log_path, std::ios::out | std::ios::trunc);
    if (!log.is_open()) {
        amrex::Print()
            << "[DumpConfigurationLog] ERROR: Cannot open "
            << log_path << "\n";
        return;
    }

    // =============================================================
    // Header
    // =============================================================
    log << "=====================================================================\n"
        << " COMPAS CONFIGURATION LOG\n"
        << " Created: " << CurrentTimeString() << "\n"
        << " Output Directory: " << plot_dir << "\n"
        << "---------------------------------------------------------------------\n"
        << " MPI Ranks: " << ParallelDescriptor::NProcs() << "\n"
        << "=====================================================================\n\n";

    // =============================================================
    // AMReX build information
    // =============================================================
    {
        log << "================ AMReX Build Information ================\n";
        log << "AMReX Version: " << amrex::Version() << "\n";

#ifdef AMREX_GIT_VERSION
        log << "AMReX Git Hash: " << AMREX_GIT_VERSION << "\n";
#endif
#ifdef AMREX_BUILD_INFO_H_
        log << "Compiler:       " << amrex::get_compiler_version() << "\n";
        log << "Build Machine:  " << amrex::buildInfoGetBuildMachine() << "\n";
        log << "Build Date:     " << amrex::buildInfoGetBuildDate() << "\n";
        log << "Compiler Flags: " << amrex::buildInfoGetBuildFlags() << "\n";
#endif
        log << "\n";
    }

    // =============================================================
    // Runtime Parameters
    // =============================================================
    log << "================ Runtime Parameters ================\n"
        << DumpRuntimeParameters() << "\n";

    // =============================================================
    // User-provided files
    // =============================================================
    for (const auto& section : sections)
    {
        const std::string& title = section.first;
        const std::string& path  = section.second;

        std::string contents = SafeLoadFileToString(path);
        if (contents.empty())
            continue;

        log << "================ " << title << " ================\n"
            << contents << "\n\n";
    }

    log << "=========================== END OF LOG ==============================\n";
    log.close();

    amrex::Print()
        << "Configuration log written to: "
        << log_path << "\n";
}

#endif // LOG_UTILS_H_
