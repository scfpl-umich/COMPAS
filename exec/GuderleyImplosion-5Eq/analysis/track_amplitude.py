# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Track the interface and the converging shock of the Guderley implosion.

For every plotfile, on the quarter plane 0 <= theta <= pi/2:

  R_eq(t)   equivalent radius of the inner gas, from its area A = int alpha_2 dA = pi R_eq^2 / 4
  h(t)      half the peak-to-valley extent of the alpha_1 = 0.5 contour, (r_max - r_min) / 2
  a_n(t)    Fourier amplitude of mode n of the interface r(theta),
                a_n = (2/Theta) int r(theta) cos(n theta) dtheta,   Theta = pi/2,
            reported while r(theta) is single-valued (before the interface rolls up)
  r_shock   smallest r on the x axis where the pressure exceeds 1.5, until the shock
            reaches the axis

The results are written to amplitude.csv and plotted.

    python analysis/track_amplitude.py plot/Guderley --mode 8 --a0 0.02
"""

import argparse
import csv

import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compas_io import find_plotfiles, load_uniform, interface_contours

parser = argparse.ArgumentParser()
parser.add_argument("plot_dir", help="directory with the plotfiles, e.g. plot/Guderley")
parser.add_argument("--mode", type=int, default=8, help="perturbation mode number n")
parser.add_argument("--a0", type=float, default=0.02, help="initial amplitude, for scaling")
parser.add_argument("--csv", default="amplitude.csv")
parser.add_argument("--out", default="amplitude.png")
args = parser.parse_args()

theta = np.linspace(0.0, 0.5 * np.pi, 721)
Theta = 0.5 * np.pi
rows = []

for plotfile in find_plotfiles(args.plot_dir):
    t, x, y, d = load_uniform(plotfile)
    dA = (x[1] - x[0]) * (y[1] - y[0])

    # Equivalent radius of the inner gas (phase 2) from its area on the quarter plane
    R_eq = np.sqrt(4.0 * np.sum(1.0 - d["a1"]) * dA / np.pi)

    # Interface: all alpha_1 = 0.5 contours, in polar coordinates
    lines = interface_contours(x, y, d["a1"])
    r_all = np.concatenate([np.hypot(l[:, 0], l[:, 1]) for l in lines])
    h = 0.5 * (r_all.max() - r_all.min())

    # Mode amplitude from the main contour, while r(theta) is single-valued
    main = lines[0]
    th = np.arctan2(main[:, 1], main[:, 0])
    single_valued = np.all(np.diff(th) > 0) or np.all(np.diff(th) < 0)
    if single_valued:
        order = np.argsort(th)
        r_theta = np.interp(theta, th[order], np.hypot(main[order, 0], main[order, 1]))
        a_n = 2.0 * trapezoid(r_theta * np.cos(args.mode * theta), theta) / Theta
    else:
        a_n = np.nan

    # Converging shock: innermost point on the x axis with P > 1.5, until it reaches the axis
    shocked = np.nonzero(d["P"][:, 0] > 1.5)[0]
    r_shock = x[shocked[0]] if shocked.size and shocked[0] > 0 else np.nan

    rows.append((t, R_eq, h, a_n, r_shock))
    print(f"t = {t:5.3f}   R_eq = {R_eq:.4f}   h = {h:.4f}   a_{args.mode} = {a_n:+.5f}   r_shock = {r_shock:.4f}")

with open(args.csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["t", "R_eq", "h", f"a_{args.mode}", "r_shock"])
    writer.writerows(rows)
print(f"Wrote {args.csv}")

t, R_eq, h, a_n, r_shock = np.array(rows).T
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)

ax1.plot(t, R_eq, "-", color="tab:blue", label=r"interface, $R_{eq}$")
ax1.plot(t, r_shock, "--", color="k", label="converging shock")
ax1.set_xlabel("t")
ax1.set_ylabel("radius")
ax1.set_ylim(bottom=0)
ax1.legend(frameon=False)

ax2.axhline(0.0, color="0.7", lw=0.8)
ax2.plot(t, h / args.a0, "-", color="tab:blue", label=r"$h/a_0$, half peak-to-valley")
ax2.plot(t, a_n / args.a0, "-", color="tab:red", label=rf"$a_{{{args.mode}}}/a_0$, mode amplitude")
ax2.set_xlabel("t")
ax2.set_ylabel("amplitude / $a_0$")
ax2.legend(frameon=False)

fig.savefig(args.out, dpi=200)
print(f"Wrote {args.out}")
