# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Interface positions, amplitudes and layer thickness of the three-material case.

Reads the interfaces.csv written in situ by UserOutputFunction (run.user_output_int = 1) for
one or more runs, and plots
  left:   x_12 and x_23 measured from the initial position of the first interface (paper, fig. 6)
  middle: the amplitudes |a_12| and |a_23| (fig. 6)
  right:  the thickness N_I of material 2 in cells (fig. 7)

    python analysis/plot_interfaces.py plot/ThreeMaterialRM/interfaces.csv
    python analysis/plot_interfaces.py L3/plot/ThreeMaterialRM/interfaces.csv L5/plot/ThreeMaterialRM/interfaces.csv \
           --labels "3 levels" "5 levels" --out interfaces.png
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("csv", nargs="+", help="interfaces.csv files")
parser.add_argument("--labels", nargs="+", help="one label per file")
parser.add_argument("--out", default="interfaces.png")
args = parser.parse_args()
labels = args.labels or [f"run {n + 1}" for n in range(len(args.csv))]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4), constrained_layout=True)
styles = ["-", "--", ":", "-."]
for n, (path, label) in enumerate(zip(args.csv, labels)):
    d = np.genfromtxt(path, delimiter=",", names=True)
    t = d["t"] * 1e6                                   # microseconds
    x0 = d["x_12"][0]
    ls = styles[n % len(styles)]
    ax1.plot(t, (d["x_12"] - x0) * 1e3, ls, color="tab:blue", label=rf"$x_{{12}}$, {label}")
    ax1.plot(t, (d["x_23"] - x0) * 1e3, ls, color="tab:red", label=rf"$x_{{23}}$, {label}")
    ax2.plot(t, np.abs(d["a_12"]) * 1e3, ls, color="tab:blue", label=rf"$|a_{{12}}|$, {label}")
    ax2.plot(t, np.abs(d["a_23"]) * 1e3, ls, color="tab:red", label=rf"$|a_{{23}}|$, {label}")
    ax3.plot(t, d["N_I"], ls, color="k", label=label)

ax1.set_ylabel("interface position [mm]")
ax2.set_ylabel("amplitude [mm]")
ax3.set_ylabel(r"$N_I$, material 2 [cells]")
for ax in (ax1, ax2, ax3):
    ax.set_xlabel(r"t [$\mu$s]")
    ax.set_xlim(0, 1000)
    ax.legend(frameon=False, fontsize=7)
ax2.set_ylim(bottom=0)
ax3.set_ylim(bottom=0)
fig.savefig(args.out, dpi=200)
print(f"Wrote {args.out}")
