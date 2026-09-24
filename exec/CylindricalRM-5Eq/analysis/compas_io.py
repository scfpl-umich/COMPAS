# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Helpers shared by the analysis scripts of the cylindrical Richtmyer-Meshkov case.

COMPAS writes AMReX plotfiles, which yt reads directly. These helpers load a plotfile,
resample it onto a uniform grid at the finest level, and extract the interface.
"""

from pathlib import Path

import numpy as np
import yt
import contourpy

yt.set_log_level(40)  # silence yt's info messages


def find_plotfiles(plot_dir):
    """Return the plotfiles in plot_dir, sorted by step number."""
    files = sorted(p for p in Path(plot_dir).glob("plt*") if p.is_dir() and ".old" not in p.name)
    if not files:
        raise SystemExit(f"No plotfiles found in {plot_dir}")
    return files


def plotfile_time(plotfile):
    """Read the simulation time from a plotfile header without loading the data."""
    lines = (Path(plotfile) / "Header").read_text().splitlines()
    nvars = int(lines[1])
    return float(lines[nvars + 3])


def load_uniform(plotfile, fields=("a1", "rho", "P")):
    """Load a plotfile and return (time, x, y, {field: 2D array}) on the finest level.

    Coarser levels are interpolated onto the finest grid, so the arrays cover the
    whole domain. Arrays are indexed [i, j], with i along x and j along y.
    """
    ds = yt.load(str(plotfile))
    level = ds.max_level
    dims = ds.domain_dimensions * ds.refine_by**level
    grid = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    data = {f: grid[("boxlib", f)][:, :, 0].d for f in fields}
    lo = ds.domain_left_edge.d
    hi = ds.domain_right_edge.d
    x = np.linspace(lo[0], hi[0], dims[0], endpoint=False) + 0.5 * (hi[0] - lo[0]) / dims[0]
    y = np.linspace(lo[1], hi[1], dims[1], endpoint=False) + 0.5 * (hi[1] - lo[1]) / dims[1]
    return float(ds.current_time), x, y, data


def interface_contours(x, y, alpha, level=0.5):
    """Return every alpha = level contour as a list of (N, 2) arrays of (x, y) points,
    longest first."""
    gen = contourpy.contour_generator(x, y, alpha.T)
    return sorted(gen.lines(level), key=len, reverse=True)


def mirror_quarter(field):
    """Mirror a quarter-plane field (x, y >= 0) about both axes onto the full plane."""
    top = np.concatenate([field[::-1, :], field], axis=0)
    return np.concatenate([top[:, ::-1], top], axis=1)
