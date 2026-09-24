#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""
===============================================================================
CSV Plotting Utility: Self-Bootstrapping, Multi-CSV, Reproducible Visualization
===============================================================================

OVERVIEW
--------
This script is a single-file plotting utility for visualizing one or more CSV
datasets with shared column structure. It is designed to be:

• Zero-setup (auto-creates a minimal virtual environment if needed)
• Deterministic and reproducible across machines
• Robust to noisy CSV data
• Suitable for production analysis and publication-quality figures

The script supports plotting any number of CSV files simultaneously, enforcing
identical headers and consistent visual semantics across datasets.

-------------------------------------------------------------------------------
EXECUTION MODEL
-------------------------------------------------------------------------------

1. Virtual Environment Bootstrapping
-----------------------------------
When executed, the script determines whether it is already running inside a
Python virtual environment.

• If NOT inside a virtual environment:
    - Creates or reuses a minimal venv at a fixed location
    - Ensures numpy and matplotlib are installed
    - Re-invokes itself inside that environment automatically

• If ALREADY inside a virtual environment:
    - Executes plotting logic directly

This guarantees dependency availability without polluting system Python.

-------------------------------------------------------------------------------
CSV INPUT RULES
-------------------------------------------------------------------------------

• One or more CSV files may be supplied via --csv
• All CSV files MUST have identical headers (names and order)
• Headers are normalized by:
    - Stripping UTF-8 BOM markers
    - Trimming surrounding whitespace
• Data is read using csv.DictReader

Row handling:
-------------
• Each row is processed independently
• Rows with non-numeric X or Y values are silently skipped
• No row-level error aborts execution

-------------------------------------------------------------------------------
PLOTTING SEMANTICS
-------------------------------------------------------------------------------

Color:
------
• Each CSV is assigned a unique color (matplotlib default cycle)
• The same CSV uses the same color across all Y variables

Line style:
-----------
• Each Y variable is assigned a unique line style
• Line styles are cycled if the number of Y variables exceeds the style list

Legend:
-------
• Single CSV:
      Y
• Multiple CSVs:
      Y (filename.csv)

-------------------------------------------------------------------------------
DIFFERENCE PLOTS (--diff)
-------------------------------------------------------------------------------

When --diff is enabled and more than one CSV is supplied:

• The first CSV is treated as the reference
• For each additional CSV:
      ΔY = CSVᵢ − CSV₀
• Difference curves:
      - Use the same color as the CSV
      - Use dot-dash ("-.") line style
      - Are labeled clearly in the legend

No interpolation is performed; datasets are truncated to common length.

-------------------------------------------------------------------------------
AXIS LIMIT OVERRIDES
-------------------------------------------------------------------------------

X-axis:
-------
• --xmin / --xmax restrict the visible X domain
• When applied, Y-axis autoscaling is recomputed using only visible data

Y-axis:
-------
• --ymin / --ymax explicitly override Y-axis limits
• If either is provided, it takes precedence over automatic Y autoscaling
• Can be combined with --xmin / --xmax

-------------------------------------------------------------------------------
COMMAND-LINE INTERFACE
-------------------------------------------------------------------------------

Required:
---------
--csv <file> [<file> ...]
--y   <col>  [<col>  ...]

Optional:
---------
--x <column>         (default: Time)
--style line|scatter (default: line)
--logx               Logarithmic X axis
--logy               Logarithmic Y axis
--xmin <value>       Minimum X-axis value
--xmax <value>       Maximum X-axis value
--ymin <value>       Minimum Y-axis value (override)
--ymax <value>       Maximum Y-axis value (override)
--diff               Plot CSVᵢ − CSV₀ differences
--list-cols          Print column names and exit

-------------------------------------------------------------------------------
OUTPUT
-------------------------------------------------------------------------------

• Interactive plot window is displayed
• A PNG file is written:
      <Y1>_<Y2>_vs_<X>.png
• Existing files with the same name are overwritten

===============================================================================
"""

import subprocess
import sys
import os
import venv
from pathlib import Path
import csv
import argparse


# ---------------------------------------------------------------------------
# Minimal virtual environment handling
# ---------------------------------------------------------------------------

def ensure_minimal_env(env_dir=None):
    """
    Ensure a minimal Python virtual environment exists with numpy and matplotlib.
    If not present, create it and install required packages.
    """
    if env_dir is None:
        script_dir = Path(__file__).resolve().parent
        env_dir = script_dir / "./compas_python_env"

    env_path = Path(env_dir).resolve()

    if not env_path.exists():
        print(f"[COMPAS] Creating minimal Python environment at '{env_path}'")
        venv.EnvBuilder(with_pip=True).create(env_path)

    if os.name == "nt":
        pip_exe = env_path / "Scripts" / "pip.exe"
        python_exe = env_path / "Scripts" / "python.exe"
    else:
        pip_exe = env_path / "bin" / "pip"
        python_exe = env_path / "bin" / "python"

    def ensure_package(pkg):
        result = subprocess.run(
            [str(python_exe), "-c", f"import {pkg}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0:
            print(f"[COMPAS] Installing {pkg}")
            subprocess.check_call([str(pip_exe), "install", pkg])

    # Keep in sync with scripts/requirements.txt
    for pkg in ("numpy", "matplotlib", "scipy", "yt", "contourpy"):
        ensure_package(pkg)

    return str(python_exe)


def running_inside_venv():
    """Return True if already running inside a virtual environment."""
    return sys.prefix != sys.base_prefix or bool(os.environ.get("VIRTUAL_ENV"))


# ---------------------------------------------------------------------------
# CSV handling
# ---------------------------------------------------------------------------

def load_csv_columns(csv_path: Path):
    """
    Load a CSV file and return normalized headers and all rows.
    """
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [
            name.strip().lstrip("\ufeff") for name in reader.fieldnames
        ]
        rows = list(reader)

    return reader.fieldnames, rows


# ---------------------------------------------------------------------------
# Plotting logic
# ---------------------------------------------------------------------------

def plot_xy(csv_paths, xvar, yvars, style, logx, logy, diff, xmin, xmax, ymin, ymax):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # --- Hard-coded styling for reproducibility ---
    mpl.rcParams.update({
        "font.family": "Times New Roman",
        "font.size": 14,
        "axes.linewidth": 1.5,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "grid.color": "#AAAAAA",
        "grid.linewidth": 0.7,
        "grid.alpha": 0.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "legend.frameon": True,
    })

    # --- Validate CSV paths ---
    for p in csv_paths:
        if not p.is_file():
            raise ValueError(f"--csv argument must be a file: {p}")

    datasets = []
    headers_ref = None

    # --- Load all CSVs and enforce identical headers ---
    for path in csv_paths:
        headers, rows = load_csv_columns(path)

        if headers_ref is None:
            headers_ref = headers
        elif headers != headers_ref:
            raise ValueError("All CSV files must have identical headers")

        datasets.append((path.name, rows))

    # --- Validate columns ---
    if xvar not in headers_ref:
        raise ValueError(f"Missing X column '{xvar}'")

    for y in yvars:
        if y not in headers_ref:
            raise ValueError(f"Missing Y column '{y}'")

    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=140)

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    line_styles = ["-", "--", ":", "-."]

    csv_colors = {
        datasets[i][0]: color_cycle[i % len(color_cycle)]
        for i in range(len(datasets))
    }

    y_linestyles = {
        yvars[i]: line_styles[i % len(line_styles)]
        for i in range(len(yvars))
    }

    # Extracted numeric data for diff computation
    extracted = {}

    # Track all plotted curves for post-xlim Y autoscaling
    plotted_curves = []

    for yi, y in enumerate(yvars):
        extracted[y] = []

        for ci, (csv_label, rows) in enumerate(datasets):
            xvals, yvals = [], []

            for r in rows:
                try:
                    xvals.append(float(r[xvar]))
                    yvals.append(float(r[y]))
                except Exception:
                    continue

            xvals = np.asarray(xvals)
            yvals = np.asarray(yvals)

            extracted[y].append((xvals, yvals))

            label = y if len(datasets) == 1 else f"{y} ({csv_label})"
            linestyle = y_linestyles[y]
            color = csv_colors[csv_label]

            if style == "scatter":
                ax.scatter(xvals, yvals, s=18, color=color, label=label)
                plotted_curves.append((xvals, yvals))
            else:
                ax.plot(
                    xvals,
                    yvals,
                    linestyle=linestyle,
                    color=color,
                    label=label,
                )
                plotted_curves.append((xvals, yvals))

        # --- Difference plots (CSV_i - CSV_0) ---
        if diff and len(datasets) >= 2:
            x_ref, y_ref = extracted[y][0]

            for di in range(1, len(extracted[y])):
                x_i, y_i = extracted[y][di]
                n = min(len(x_ref), len(x_i))

                ydiff = y_i[:n] - y_ref[:n]

                ax.plot(
                    x_ref[:n],
                    ydiff,
                    linestyle="-.",
                    color=csv_colors[datasets[di][0]],
                    label=f"{y} (Δ {datasets[di][0]})",
                )

                plotted_curves.append((x_ref[:n], ydiff))

    # --- Axis scaling ---
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    # Apply X limits and recompute Y limits from visible data only
    if xmin is not None or xmax is not None:
        ax.set_xlim(left=xmin, right=xmax)

        ymin_auto, ymax_auto = None, None

        for xvals, yvals in plotted_curves:
            mask = np.ones_like(xvals, dtype=bool)

            if xmin is not None:
                mask &= xvals >= xmin
            if xmax is not None:
                mask &= xvals <= xmax

            if not np.any(mask):
                continue

            y_visible = yvals[mask]

            if logy:
                y_visible = y_visible[y_visible > 0]

            if y_visible.size == 0:
                continue

            ymin_auto = y_visible.min() if ymin_auto is None else min(ymin_auto, y_visible.min())
            ymax_auto = y_visible.max() if ymax_auto is None else max(ymax_auto, y_visible.max())

        if ymin_auto is not None and ymax_auto is not None:
            ax.set_ylim(ymin_auto, ymax_auto)

    # Explicit Y overrides take precedence over autoscaling
    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)

    ax.set_xlabel(xvar)
    ax.set_ylabel(", ".join(yvars))
    ax.set_title(f"{', '.join(yvars)} vs {xvar}")
    ax.grid(True)

    ax.legend(loc="upper right", fontsize=12)
    fig.tight_layout()

    outname = f"{'_'.join(yvars)}_vs_{xvar}.png"
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Flexible multi-CSV plotting utility")

    parser.add_argument("--csv", type=Path, nargs="+", required=True)
    parser.add_argument("--x", default="Time")
    parser.add_argument("--y", nargs="+")
    parser.add_argument("--style", choices=["line", "scatter"], default="line")
    parser.add_argument("--logx", action="store_true")
    parser.add_argument("--logy", action="store_true")
    parser.add_argument("--xmin", type=float)
    parser.add_argument("--xmax", type=float)
    parser.add_argument("--ymin", type=float)
    parser.add_argument("--ymax", type=float)
    parser.add_argument("--diff", action="store_true")
    parser.add_argument("--list-cols", action="store_true")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    if running_inside_venv():
        headers, _ = load_csv_columns(args.csv[0])

        if args.list_cols:
            print("Available columns:")
            for h in headers:
                print(f"  - {h}")
            sys.exit(0)

        if not args.y:
            raise SystemExit("Error: --y must be specified unless --list-cols is used")

        plot_xy(
            csv_paths=args.csv,
            xvar=args.x,
            yvars=args.y,
            style=args.style,
            logx=args.logx,
            logy=args.logy,
            diff=args.diff,
            xmin=args.xmin,
            xmax=args.xmax,
            ymin=args.ymin,
            ymax=args.ymax,
        )

    else:
        python_exe = ensure_minimal_env()

        env = os.environ.copy()
        env["VIRTUAL_ENV"] = str(Path(python_exe).parent.parent)
        env["PATH"] = str(Path(python_exe).parent) + os.pathsep + env["PATH"]

        subprocess.check_call(
            [python_exe, __file__, *sys.argv[1:]],
            env=env,
        )
