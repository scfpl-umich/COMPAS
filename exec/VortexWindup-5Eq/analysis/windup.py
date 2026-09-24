# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Errors of the wound-up interface against the exact solution.

The flow is steady and both materials are the same gas, so a point at radius r from the center
turns at Omega(r) = beta/(2 pi rc) exp((1 - r^2/rc^2)/2), and the exact alpha_2 at time t is the
initial one turned by -Omega(r) t. For each plotfile of each run the script averages the exact
alpha_2 over ns x ns points per cell (on the finest grid; coarser levels are injected) and prints

  t          time
  L1         int |alpha_2 - alpha_2,exact| dA
  L1_cells   L1 / (P dx), with P the length of the exact interface (the total variation of the
             exact cell averages) and dx the finest spacing: the mean width of the error, in cells
  pieces_k   connected pieces of {alpha_k > 0.5} (8-neighbour), material k = 1, 2
  exact_k    the same for the exact alpha_k
  extra      extra pieces, max(pieces_1 - exact_1, 0) + max(pieces_2 - exact_2, 0): the beads
             and specks the arms break into
  specks     pieces of fewer than 4 cells, both materials
  mixed      cells with 0.01 < alpha_1 < 0.99
  a1_min, a1_max   extremes of alpha_1

The problem parameters are read from the inputs file (prob/inputs by default), with the same
KEY=VALUE overrides as the run. With --plot, the script also draws the exact alpha_2 and each run
at one time (--time; by default the last plotfile of each run) in a window about the vortex.

    python analysis/windup.py plot/VortexWindup
    python analysis/windup.py old/plot/VortexWindup new/plot/VortexWindup --csv windup.csv \\
           --plot windup.png --time 20
    python analysis/windup.py plot/VortexWindup --set prob.width=0.5
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np
import yt
from scipy import ndimage

yt.set_log_level(40)  # silence yt's info messages

HERE = Path(__file__).resolve().parent
STYLE = HERE.parents[2] / "scripts" / "compas.mplstyle"   # the style of the documentation figures

# Teal to sand, the map of the volume-fraction figures of the documentation: OKLCH anchors
# (L, C, h), with L rising linearly and C, h interpolated piecewise in L. Material 1 is at the
# teal end and material 2 at the sand end.
ANCHORS = np.array([(0.30, 0.085, 252), (0.50, 0.080, 215), (0.72, 0.060, 150), (0.86, 0.070, 95),
                    (0.95, 0.045, 85)])


def teal_sand(n=256):
    """sRGB table of the teal-sand map, n entries."""
    M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929], [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468], [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    L = np.linspace(ANCHORS[0, 0], ANCHORS[-1, 0], n)
    C = np.interp(L, ANCHORS[:, 0], ANCHORS[:, 1])
    h = np.radians(np.interp(L, ANCHORS[:, 0], ANCHORS[:, 2]))
    lab = np.stack([L, C * np.cos(h), C * np.sin(h)], -1)
    lin = ((lab @ np.linalg.inv(M2).T) ** 3) @ np.linalg.inv(M1).T
    return np.clip(np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * np.abs(lin) ** (1 / 2.4) - 0.055), 0, 1)


# ---------------------------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------------------------

DEFAULTS = {"beta": 5.0, "rc": 1.0, "xc": 0.0, "yc": 0.0, "angle": 0.0, "offset": 0.0, "width": 0.0}


def read_params(inputs, overrides):
    """prob.* values of the inputs file, with KEY=VALUE overrides, as floats."""
    p = dict(DEFAULTS)
    lines = Path(inputs).read_text().splitlines() + list(overrides)
    for line in lines:
        m = re.match(r"\s*prob\.(\w+)\s*=\s*([^#\s]+)", line.split("#")[0])
        if m and m.group(1) in p:
            p[m.group(1)] = float(m.group(2))
    return p


# ---------------------------------------------------------------------------------------------
# Plotfiles
# ---------------------------------------------------------------------------------------------

def plotfiles(plot_dir):
    """Plotfiles of a plot directory, in time order."""
    files = sorted(q for q in Path(plot_dir).glob("*[0-9]") if (q / "Header").is_file() and ".old" not in q.name)
    if not files:
        raise SystemExit(f"No plotfiles found in {plot_dir}")
    return files


def plotfile_time(plotfile):
    """Simulation time of a plotfile, from its header."""
    lines = (Path(plotfile) / "Header").read_text().splitlines()
    return float(lines[int(lines[1]) + 3])


def load(plotfile):
    """(t, x, y, dx, a1) on the finest grid of the whole domain; coarser levels are injected."""
    ds = yt.load(str(plotfile))
    lo, hi = ds.domain_left_edge.d, ds.domain_right_edge.d
    n = ds.domain_dimensions * ds.refine_by ** ds.max_level
    g = ds.covering_grid(level=ds.max_level, left_edge=lo, dims=[n[0], n[1], 1])
    a1 = g[("boxlib", "a1")][:, :, 0].d
    dx = (hi[0] - lo[0]) / n[0]
    x = lo[0] + (np.arange(n[0]) + 0.5) * dx
    y = lo[1] + (np.arange(n[1]) + 0.5) * dx
    return float(ds.current_time), x, y, dx, a1


# ---------------------------------------------------------------------------------------------
# Exact solution and metrics
# ---------------------------------------------------------------------------------------------

def exact_alpha2(p, x, y, dx, t, ns=8):
    """Exact alpha_2 at time t, averaged over ns x ns points of each cell, as [i, j]."""
    ph = math.radians(p["angle"])
    nx, ny = math.cos(ph), math.sin(ph)
    Xc, Yc = np.meshgrid(x - p["xc"], y - p["yc"], indexing="ij")
    A = np.zeros_like(Xc)
    for a in range(ns):
        for b in range(ns):
            X = Xc + ((a + 0.5) / ns - 0.5) * dx
            Y = Yc + ((b + 0.5) / ns - 0.5) * dx
            th = -p["beta"] / (2 * math.pi * p["rc"]) * np.exp(0.5 * (1 - (X * X + Y * Y) / p["rc"] ** 2)) * t
            d = nx * (np.cos(th) * X - np.sin(th) * Y) + ny * (np.sin(th) * X + np.cos(th) * Y) - p["offset"]
            A += (d > 0) & ((p["width"] <= 0) | (d < p["width"]))
    return A / ns ** 2


def pieces(mask):
    """Number of 8-connected pieces of a boolean field, and their sizes in cells."""
    lab, n = ndimage.label(mask, structure=np.ones((3, 3), int))
    return n, np.bincount(lab.ravel())[1:]


def metrics(t, dx, a1, e2):
    a2 = 1.0 - a1
    dA = dx * dx
    gx, gy = np.gradient(e2, dx)
    P = np.hypot(gx, gy).sum() * dA                       # length of the exact interface
    L1 = np.abs(a2 - e2).sum() * dA
    r = dict(t=round(t, 6), L1=L1, L1_cells=L1 / (P * dx) if P > 0 else 0.0)
    extra, specks = 0, 0
    for k, (a, e) in enumerate([(a1, 1.0 - e2), (a2, e2)], start=1):
        n, s = pieces(a > 0.5)
        ne, _ = pieces(e > 0.5)
        r[f"pieces_{k}"], r[f"exact_{k}"] = n, ne
        extra += max(n - ne, 0)
        specks += int((s < 4).sum())
    r["extra"], r["specks"] = extra, specks
    r["mixed"] = int(((a1 > 0.01) & (a1 < 0.99)).sum())
    r["a1_min"], r["a1_max"] = float(a1.min()), float(a1.max())
    return r


# ---------------------------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------------------------

def plot(frames, labels, out, half):
    """The exact alpha_2 and each run, side by side, in the window |x|, |y| < half about the center."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    plt.style.use(STYLE)
    cmap = ListedColormap(teal_sand())
    t, x, y, dx, _, e2, p = frames[0]
    fields = [e2] + [1.0 - f[4] for f in frames]
    n = len(fields)
    fig, axs = plt.subplots(1, n, figsize=(2.6 * n, 2.75), sharey=True, squeeze=False)
    ext = [x[0] - dx / 2, x[-1] + dx / 2, y[0] - dx / 2, y[-1] + dx / 2]
    for k, (ax, f) in enumerate(zip(axs[0], fields)):
        ax.imshow(np.clip(f, 0, 1).T, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=1,
                  interpolation="nearest")
        ax.set_xlim(p["xc"] - half, p["xc"] + half)
        ax.set_ylim(p["yc"] - half, p["yc"] + half)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$")
        ax.minorticks_off()
    axs[0][0].set_ylabel(r"$y$")
    fig.subplots_adjust(wspace=0.06)
    fig.savefig(out, dpi=300)
    print(f"Wrote {out} (panels: exact, {', '.join(labels)}; t = {t:g})")


# ---------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("plot_dirs", nargs="+", help="plot directories, one per run")
    ap.add_argument("--labels", nargs="+", help="one label per run (default: the directory names)")
    ap.add_argument("--inputs", default=str(HERE.parent / "prob" / "inputs"), help="inputs file of the runs")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE", help="overrides of the runs, e.g. prob.width=0.5")
    ap.add_argument("--ns", type=int, default=8, help="samples per direction of the exact solution in each cell")
    ap.add_argument("--csv", help="write the metrics to this file")
    ap.add_argument("--plot", help="write a figure of the exact and computed alpha_2 to this file")
    ap.add_argument("--time", type=float, help="time of the figure (default: the last plotfile of the first run)")
    ap.add_argument("--window", type=float, default=3.0, help="half width of the window of the figure")
    a = ap.parse_args()

    p = read_params(a.inputs, a.set)
    labels = a.labels or list(a.plot_dirs)
    rows, frames = [], []
    for d, label in zip(a.plot_dirs, labels):
        files = plotfiles(d)
        pick = None
        if a.plot:
            if a.time is None:
                pick = files[-1]
            else:
                times = [plotfile_time(q) for q in files]
                pick = files[int(np.argmin(np.abs(np.array(times) - a.time)))]
        for q in files:
            t, x, y, dx, a1 = load(q)
            e2 = exact_alpha2(p, x, y, dx, t, a.ns)
            r = dict(run=label, plotfile=q.name)
            r.update(metrics(t, dx, a1, e2))
            rows.append(r)
            print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()), flush=True)
            if q == pick:
                frames.append((t, x, y, dx, a1, e2, p))
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {a.csv}")
    if a.plot:
        plot(frames, labels, a.plot, a.window)


if __name__ == "__main__":
    sys.exit(main())
