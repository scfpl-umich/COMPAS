#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""
COMPAS convergence utility

Run from the directory of a case that defines CONVERGENCE and exact_solution. At the end of
each run COMPAS writes the error file <conv_output_file>.0 (see COMPAS::OutputError):

    <conv_prob_name> Nx Ny dt runtime
    <variable> L2 Linf          (one line per conservative variable)

conv_prob_name is empty unless convergence.conv_prob_name is set.

--run            run the case at each resolution and save the errors of every variable
--plot           plot the saved errors on log-log axes, with a line through the last two points,
                 and save the plot next to the CSV as a PNG with the same name
--res N ...      resolutions; each run sets amr.n_cell=N N 1 (default 8 16 32 64 128 256)
--norm NORM      Linf (default) or L2, the norm saved by --run and named on the plot
--exe EXE        executable (default: the only *.ex in the current directory)
--inputs FILE    inputs file (default ./prob/inputs)
--nprocs N       MPI ranks (default 1); with 1 the executable runs without mpirun
--csv FILE       file written by --run and read by --plot (default convergence_all.csv)
--vars VAR ...   variables to plot (default all); --run always saves every variable
--clear BOOL     true (default) runs exec/clean.sh in the case directory before the runs

The CSV has the columns h, var1, var2, ..., where h = (prob_hi - prob_lo)/N is the grid
spacing in x, with the domain read from geometry.prob_lo and geometry.prob_hi in the inputs
file. The order between successive resolutions is printed after --run.
"""

import subprocess
import sys
import os
import venv
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CSV = "convergence_all.csv"


# ---------------------------------------------------------------------------
# Minimal virtual environment handling
# ---------------------------------------------------------------------------

def ensure_minimal_env(env_dir=None):
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
            subprocess.check_call([str(pip_exe), "install", pkg])

    # Keep in sync with scripts/requirements.txt
    for pkg in ("numpy", "matplotlib", "scipy", "yt", "contourpy"):
        ensure_package(pkg)

    return str(python_exe)


def running_inside_venv():
    return sys.prefix != sys.base_prefix or bool(os.environ.get("VIRTUAL_ENV"))


# ---------------------------------------------------------------------------
# Inputs file and executable
# ---------------------------------------------------------------------------

def read_domain_length(inputs):
    """Length of the domain in x, from the last geometry.prob_lo and prob_hi in the inputs."""
    values = {}
    with open(inputs) as f:
        for line in f:
            line = line.split("#", 1)[0]
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            key = key.strip()
            if key in ("geometry.prob_lo", "geometry.prob_hi"):
                values[key] = float(val.split()[0])
    if len(values) != 2:
        sys.exit(f"[ERROR] geometry.prob_lo and geometry.prob_hi not found in {inputs}")
    return values["geometry.prob_hi"] - values["geometry.prob_lo"]


def find_executable():
    exes = sorted(p.name for p in Path(".").glob("*.ex")
                  if p.is_file() and os.access(p, os.X_OK))
    if len(exes) != 1:
        sys.exit(f"[ERROR] Found {len(exes)} executables (*.ex) here; choose one with --exe")
    return "./" + exes[0]


# ---------------------------------------------------------------------------
# Run Convergence (ALL VARIABLES)
# ---------------------------------------------------------------------------

def run_convergence(
    exe,
    base_inputs,
    resolutions,
    nprocs=1,
    conv_file="conv_temp",
    save_csv=DEFAULT_CSV,
    clear="true",
    norm="Linf"
):
    import numpy as np
    import csv

    hs = []
    all_errors = {}
    variable_names = None
    column = {"L2": 1, "Linf": 2}[norm]
    length = read_domain_length(base_inputs)

    if clear == "true":
        subprocess.run(["bash", str(SCRIPT_DIR.parent / "exec" / "clean.sh")], check=True)

    for N in resolutions:
        print(f"\nRunning resolution {N}x{N} ...")

        # COMPAS writes the errors to conv_file + ".0"; remove it so a run that does not
        # write it cannot leave the errors of the previous resolution in its place
        err_file = conv_file + ".0"
        if os.path.exists(err_file):
            os.remove(err_file)

        run_cmd = [
            exe,
            base_inputs,
            f"amr.n_cell={N} {N} 1",
            f"convergence.conv_output_file={conv_file}"
        ]
        if nprocs > 1:
            run_cmd = ["mpirun", "-n", str(nprocs)] + run_cmd
        subprocess.run(run_cmd, check=True)

        if not os.path.exists(err_file):
            sys.exit(f"[ERROR] {err_file} was not written; is CONVERGENCE defined for this case?")
        with open(err_file, "r") as f:
            lines = f.readlines()

        # First line: <conv_prob_name> Nx Ny dt runtime, where the name may be empty
        Nx = int(lines[0].split()[-4])
        if Nx != N:
            sys.exit(f"[ERROR] {err_file} is for Nx = {Nx}, not {N}")
        h = length / Nx
        hs.append(h)

        current_vars = []
        current_errors = []

        for line in lines[1:]:
            tokens = line.strip().split()
            if not tokens:
                continue
            current_vars.append(tokens[0])
            current_errors.append(float(tokens[column]))

        if variable_names is None:
            variable_names = current_vars
            for v in variable_names:
                all_errors[v] = []

        for v, e in zip(current_vars, current_errors):
            all_errors[v].append(e)

    hs = np.array(hs)

    # Save CSV
    with open(save_csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["h"] + variable_names
        writer.writerow(header)
        for i in range(len(hs)):
            row = [hs[i]] + [all_errors[v][i] for v in variable_names]
            writer.writerow(row)

    print(f"\n[INFO] Saved the {norm} errors of all variables to {save_csv}")

    # ---------------------------------------------------------
    # Compute and print the order between successive resolutions
    # ---------------------------------------------------------
    print(f"\n{norm} order between successive resolutions (the last is the final rate):")
    print("-------------------------------------------------")
    print(f"{'N':20s} : " + " ".join(f"{a:>4d}-{b:<4d}" for a, b in
                                      zip(resolutions[:-1], resolutions[1:])))

    for var in variable_names:
        errors = np.array(all_errors[var])

        if len(errors) < 2:
            print(f"{var:20s} : ERROR computing rate (Not enough data points)")
            continue

        rates = []
        for i in range(1, len(errors)):
            if errors[i] <= 0 or errors[i-1] <= 0:
                rates.append(f"{'n/a':>9s}")
                continue
            slope = (np.log(errors[i]) - np.log(errors[i-1])) / \
                    (np.log(hs[i]) - np.log(hs[i-1]))
            rates.append(f"{slope:9.4f}" if np.isfinite(slope) else f"{'n/a':>9s}")

        print(f"{var:20s} : " + " ".join(rates))

    print("-------------------------------------------------\n")

    return hs, all_errors


# ---------------------------------------------------------------------------
# Plot (ALL or SUBSET)
# ---------------------------------------------------------------------------

def plot_convergence(csv_file, selected_vars=None, norm="Linf"):
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.genfromtxt(csv_file, delimiter=",", names=True)
    hs = data["h"]

    all_vars = list(data.dtype.names)
    all_vars.remove("h")

    if selected_vars:
        vars_to_plot = selected_vars
    else:
        vars_to_plot = all_vars

    plt.figure(figsize=(7,6))

    for var in vars_to_plot:
        errors = data[var]

        plt.loglog(hs, errors, "o-", label=var)

        try:
            if len(errors) < 2:
                raise ValueError("Not enough data points")

            if errors[-1] <= 0 or errors[-2] <= 0:
                raise ValueError("Non-positive error encountered")

            slope = (np.log(errors[-1]) - np.log(errors[-2])) / \
                    (np.log(hs[-1]) - np.log(hs[-2]))

            if not np.isfinite(slope):
                raise ValueError("Slope is not finite")

            # Fit line using last two points
            h_fit = np.array([hs[-2], hs[-1]])
            e_fit = errors[-1] * (h_fit / hs[-1])**slope

            plt.loglog(
                h_fit,
                e_fit,
                "--",
                linewidth=1.5,
                label=f"{var} fit (p≈{slope:.2f})"
            )

        except Exception as e:
            print(f"[WARNING] Could not compute slope for {var}: {e}")

    plt.gca().invert_xaxis()
    plt.xlabel("Grid spacing h")
    plt.ylabel(f"{norm} error")
    plt.title("Convergence Rates")
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    png = Path(csv_file).with_suffix(".png")
    plt.savefig(png, dpi=200)
    print(f"[INFO] Saved the plot to {png}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Run a COMPAS case over several resolutions and plot its order of accuracy")
    parser.add_argument("--exe",
                        help="Executable (default: the only *.ex in the current directory)")
    parser.add_argument("--inputs", default="./prob/inputs", help="Inputs file")
    parser.add_argument("--res", type=int, nargs="+",
                        default=[8,16,32,64,128,256],
                        help="Resolutions; each run sets amr.n_cell=N N 1")
    parser.add_argument("--nprocs", type=int, default=1,
                        help="MPI ranks; with 1 the executable runs without mpirun")
    parser.add_argument("--norm", choices=["Linf", "L2"], default="Linf",
                        help="Error norm saved by --run and named on the plot")
    parser.add_argument("--csv", default=DEFAULT_CSV,
                        help="File written by --run and read by --plot")
    parser.add_argument("--vars", nargs="+",
                        help="Subset of variables to plot (--run saves every variable)")
    parser.add_argument("--run", action="store_true",
                        help="Run the case at each resolution and save the errors")
    parser.add_argument("--plot", action="store_true",
                        help="Plot the saved errors")
    parser.add_argument("--clear", choices=["true", "false"], default="true",
                        help="Run exec/clean.sh in the case directory before the runs")

    args = parser.parse_args()
    if not (args.run or args.plot):
        parser.error("nothing to do; give --run, --plot or both")
    return args


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    if running_inside_venv():

        if args.run:
            run_convergence(
                exe=args.exe or find_executable(),
                base_inputs=args.inputs,
                resolutions=args.res,
                nprocs=args.nprocs,
                save_csv=args.csv,
                clear=args.clear,
                norm=args.norm
            )

        if args.plot:
            plot_convergence(args.csv, args.vars, args.norm)

    else:
        python_exe = ensure_minimal_env()
        env = os.environ.copy()
        env["VIRTUAL_ENV"] = str(Path(python_exe).parent.parent)
        env["PATH"] = str(Path(python_exe).parent) + os.pathsep + env["PATH"]

        subprocess.check_call(
            [python_exe, __file__, *sys.argv[1:]],
            env=env,
        )