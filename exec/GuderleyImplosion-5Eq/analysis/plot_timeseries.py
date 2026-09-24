# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Plot the in-situ interface history written by UserOutputFunction.

interface.csv has one line per coarse step (run.user_output_int = 1). The shock is compared
with Guderley's R = Rs (1 - t/tau_c)^alpha from similarity.csv, which it follows until it
reaches the interface. If the plotfile-based results of track_amplitude.py are given with
--plotfiles, they are drawn as markers on top, which checks the two against each other.

    python analysis/plot_timeseries.py plot/Guderley/interface.csv --plotfiles amplitude.csv
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("csv", help="interface.csv written by UserOutputFunction")
parser.add_argument("--plotfiles", help="amplitude.csv from track_amplitude.py (optional)")
parser.add_argument("--mode", type=int, default=8)
parser.add_argument("--a0", type=float, default=0.02)
parser.add_argument("--out", default="timeseries.png")
args = parser.parse_args()

d = np.genfromtxt(args.csv, delimiter=",", names=True)
a_n = d[f"a_{args.mode}"]

# The shock radius is meaningful until the shock reaches the origin
focus = np.argmax(d["r_shock"] < 0.01)
r_shock = np.where(np.arange(d.size) < focus, d["r_shock"], np.nan)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)

ax1.plot(d["t"], d["R_eq"], "-", color="tab:blue", label=r"interface, $R_{eq}$")
ax1.plot(d["t"], r_shock, "--", color="k", label="converging shock")

# Guderley's trajectory without the interface
ss = np.genfromtxt(Path(args.csv).parent / "similarity.csv", delimiter=",", names=True)
tg = np.linspace(0.0, float(ss["tau_c"]), 200)
ax1.plot(tg, ss["Rs"] * (1.0 - tg / ss["tau_c"])**ss["alpha"], ":", color="0.4",
         label=rf"Guderley, $\alpha$ = {float(ss['alpha']):.4f}")
ax1.set_xlabel("t")
ax1.set_ylabel("radius")
ax1.set_ylim(bottom=0)

ax2.axhline(0.0, color="0.7", lw=0.8)
ax2.plot(d["t"], a_n / args.a0, "-", color="tab:red", label=rf"$a_{{{args.mode}}}/a_0$, mode amplitude")
ax2.plot(d["t"], d["h"] / args.a0, "-", color="tab:blue", label=r"$h/a_0$, half peak-to-valley")
ax2.set_xlabel("t")
ax2.set_ylabel("amplitude / $a_0$")

if args.plotfiles:
    p = np.genfromtxt(args.plotfiles, delimiter=",", names=True)
    ax1.plot(p["t"], p["R_eq"], "o", ms=3, mfc="none", color="tab:blue", label="from plotfiles")
    ax2.plot(p["t"], p[f"a_{args.mode}"] / args.a0, "o", ms=3, mfc="none", color="tab:red")
    ax2.plot(p["t"], p["h"] / args.a0, "o", ms=3, mfc="none", color="tab:blue", label="from plotfiles")

ax1.legend(frameon=False, fontsize=8)
ax2.legend(frameon=False, fontsize=8)
fig.savefig(args.out, dpi=200)
print(f"Wrote {args.out}  ({d.size} samples)")
