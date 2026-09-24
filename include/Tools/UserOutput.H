// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef USER_OUTPUT_H_
#define USER_OUTPUT_H_

// =======================================
// UserOutput: in-situ diagnostics for a case
//
// A case that defines USER_OUTPUT_FUNC in prob/ProblemICBC.H provides
//     template <class T> void UserOutputFunction(T& o)
// which COMPAS calls at the start of the run and then every run.user_output_int coarse steps,
// or with each plotfile when run.user_output_int <= 0. These helpers let the function state
// only what it measures. Every rank calls them: Measure is collective, and WriteRow writes
// from the I/O rank.
//
// Measure(o, ops, ngrow, f) evaluates the GPU lambda f on the composite grid, the cells of
// each level that no finer level covers, and combines each of its N outputs over the cells
// and the MPI ranks with the matching operation in ops: Sum, Min or Max. Measure<N>(o, op,
// ngrow, f) applies one operation to all N outputs, and a last argument FinestLevel evaluates
// f on the finest level only. For each cell, f(c, U, parm) returns Values<N> from
//   c.i, c.j, c.k   the index of the cell on its level
//   c.x, c.dx       its center and size
//   c.dV            its volume, an area in 2D
//   c.domain        the index box of the domain on its level
//   c.finest        true on the finest level
//   U               the conservative state U(i,j,k,n), filled up to ngrow cells around c
//   parm            the case's Parm
// A Sum adds the values themselves, so an integral returns its integrand times c.dV. A cell
// with nothing to contribute returns 0 for a Sum, a large value such as 1e30 for a Min, and a
// large negative value such as -1e30 for a Max.
//
// WriteRow(o, file, columns) appends the time, the step and the named values as one line of
// <output directory>/<case name>/<file>, and at step 0 starts the file with a header. Step(o)
// and Time(o) are the coarse step and the time, HostParm(o) is the case's Parm on the host,
// FinestGeom(o) the geometry of the finest level and Directory(o) the directory of the files.
//
// Example: the volume of material 2, and the largest x of a cell where it is the majority
//
//   #define USER_OUTPUT_FUNC
//   template <class T>
//   void UserOutputFunction(T& o)
//   {
//       using namespace UserOutput;
//       const auto v = Measure(o, {Sum, Max}, 0,
//       [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<2>
//       {
//           const Real a2 = 1.0 - U(c.i, c.j, c.k, INDEX_VolumeFraction1);
//           return {a2 * c.dV, (a2 >= 0.5) ? c.x[0] : -1.0e30};
//       });
//       WriteRow(o, "material2.csv", {{"V_2", v[0]}, {"x_2", v[1]}});
//   }
// =======================================

#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

namespace UserOutput
{
using amrex::Real;

enum Op     { Sum, Min, Max };
enum Levels { AllLevels, FinestLevel };

using State = amrex::Array4<Real const>;
template <int N> using Values = amrex::GpuArray<Real, N>;

struct Cell
{
    int i, j, k;
    amrex::GpuArray<Real, AMREX_SPACEDIM> x, dx;   // center and size
    Real       dV;                                 // volume, an area in 2D
    amrex::Box domain;                             // index box of the domain on this level
    bool       finest;                             // true on the finest level
};

template <class T> int         Step       (T const& o) { return o.istep[0]; }
template <class T> Real        Time       (T const& o) { return o.t_new[0]; }
template <class T> auto const& HostParm   (T const& o) { return *o.h_parm; }
template <class T> auto const& FinestGeom (T const& o) { return o.Geom(o.finestLevel()); }
template <class T> std::string Directory  (T const& o) { return o.output_dir + "/" + o.case_name + "/"; }


template <int N, class T, class F>
Values<N> Measure (T& o, Op const (&ops)[N], int ngrow, F const& f, Levels levels = AllLevels)
{
    const Real hi = std::numeric_limits<Real>::max(), lo = std::numeric_limits<Real>::lowest();
    amrex::GpuArray<int, N> op;
    Values<N> r_sum, r_min, r_max;
    for (int n = 0; n < N; ++n) { op[n] = ops[n]; r_sum[n] = 0.0; r_min[n] = hi; r_max[n] = lo; }

    auto const* parm   = o.d_parm;
    const int   finest = o.finestLevel();
    for (int lev = (levels == FinestLevel) ? finest : 0; lev <= finest; ++lev)
    {
        // The state, with ngrow ghost cells filled from this level, the coarser ones and the boundaries
        amrex::MultiFab const& S = o.dof_new[lev];
        amrex::MultiFab G;
        if (ngrow > 0) {
            G.define(S.boxArray(), S.DistributionMap(), S.nComp(), ngrow);
            o.FillPatch(G, lev, o.t_new[0], o.dof_new, o.t_new, o.dof_new, o.t_new, 0, 0, S.nComp());
        }
        amrex::MultiFab const& U = (ngrow > 0) ? G : S;

        // 1 on the cells that the next finer level does not cover
        amrex::iMultiFab mask = (lev < finest)
            ? amrex::makeFineMask(S.boxArray(), S.DistributionMap(), o.dof_new[lev+1].boxArray(), o.refRatio(lev), 1, 0)
            : amrex::iMultiFab(S.boxArray(), S.DistributionMap(), 1, 0);
        if (lev == finest) { mask.setVal(1); }

        // The outputs of f on the uncovered cells, and the identity of each operation on the covered ones
        const auto       dx        = o.Geom(lev).CellSizeArray();
        const auto       plo       = o.Geom(lev).ProbLoArray();
        const amrex::Box domain    = o.Geom(lev).Domain();
        const Real       dV        = AMREX_D_TERM(dx[0], * dx[1], * dx[2]);
        const bool       on_finest = (lev == finest);
        amrex::MultiFab Q(S.boxArray(), S.DistributionMap(), N, 0);
        for (amrex::MFIter mfi(Q, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            auto const& u = U.const_array(mfi);
            auto const& m = mask.const_array(mfi);
            auto const& q = Q.array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Values<N> v;
                if (m(i,j,k)) {
                    const Cell c{i, j, k, {AMREX_D_DECL(plo[0] + (i + 0.5) * dx[0],
                                                        plo[1] + (j + 0.5) * dx[1],
                                                        plo[2] + (k + 0.5) * dx[2])},
                                 dx, dV, domain, on_finest};
                    v = f(c, u, *parm);
                } else {
                    for (int n = 0; n < N; ++n) { v[n] = (op[n] == Sum) ? Real(0.0) : ((op[n] == Min) ? hi : lo); }
                }
                for (int n = 0; n < N; ++n) { q(i,j,k,n) = v[n]; }
            });
        }
        for (int n = 0; n < N; ++n) {
            if      (op[n] == Sum) { r_sum[n] += Q.sum(n, true); }
            else if (op[n] == Min) { r_min[n] = amrex::min(r_min[n], Q.min(n, 0, true)); }
            else                   { r_max[n] = amrex::max(r_max[n], Q.max(n, 0, true)); }
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(r_sum.data(), N);
    amrex::ParallelDescriptor::ReduceRealMin(r_min.data(), N);
    amrex::ParallelDescriptor::ReduceRealMax(r_max.data(), N);
    Values<N> r;
    for (int n = 0; n < N; ++n) { r[n] = (op[n] == Sum) ? r_sum[n] : ((op[n] == Min) ? r_min[n] : r_max[n]); }
    return r;
}

template <int N, class T, class F>
Values<N> Measure (T& o, Op op, int ngrow, F const& f, Levels levels = AllLevels)
{
    Op ops[N];
    for (int n = 0; n < N; ++n) { ops[n] = op; }
    return Measure(o, ops, ngrow, f, levels);
}


template <class T>
void WriteRow (T const& o, std::string const& file, std::vector<std::pair<std::string, Real>> const& columns)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    const std::string path = Directory(o) + file;
    std::ofstream out(path, (Step(o) == 0) ? std::ios::trunc : std::ios::app);
    if (!out) { amrex::Warning("UserOutput::WriteRow cannot open " + path); return; }
    out.precision(std::numeric_limits<Real>::max_digits10);
    if (Step(o) == 0) {
        out << "t,step";
        for (auto const& col : columns) { out << "," << col.first; }
        out << "\n";
    }
    out << Time(o) << "," << Step(o);
    for (auto const& col : columns) { out << "," << col.second; }
    out << "\n";
}

}  // namespace UserOutput

#endif
