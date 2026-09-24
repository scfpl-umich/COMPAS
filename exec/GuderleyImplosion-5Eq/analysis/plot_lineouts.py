# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Density, velocity and pressure along the x axis against Guderley's similarity solution.

Run the case without the interface (prob.R0=0 prob.a0=0), so that the shock converges on
the axis as in Guderley's problem, and compare the plotfiles before the collapse time
tau_c with the exact solution tabulated by the code (similarity*.csv, written at step 0):
    u = -(r/(tau_c - t)) U(xi),  c^2 = (r/(tau_c - t))^2 C(xi),  rho = rho0 G(xi),
    xi = r/R(t),  R(t) = Rs (1 - t/tau_c)^alpha,  and the w-table far behind the shock.

    mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs prob.R0=0 prob.a0=0 stop_time=0.25 amr.case_name=./plot/NoInterface
    python analysis/plot_lineouts.py plot/NoInterface --times 0 0.08 0.16 0.24
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

from compas_io import find_plotfiles, plotfile_time

yt.set_log_level(40)

parser = argparse.ArgumentParser()
parser.add_argument("plot_dir")
parser.add_argument("--times", type=float, nargs="+", default=[0.0, 0.08, 0.16, 0.24])
parser.add_argument("--out", default="lineouts.png")
args = parser.parse_args()

d = Path(args.plot_dir)
ss = np.genfromtxt(d / "similarity.csv", delimiter=",", names=True)
tx = np.genfromtxt(d / "similarity_xi.csv", delimiter=",", names=True)
tw = np.genfromtxt(d / "similarity_w.csv", delimiter=",", names=True)
a, tau_c, Rs, A1a, g, rho0, P0 = (float(ss[k]) for k in ("alpha", "tau_c", "Rs", "A1a", "gamma", "rho0", "P0"))


def exact(r, t):
    """rho, u_r and P of the similarity solution at radii r and time t < tau_c."""
    dt = tau_c - t
    R = Rs * (dt / tau_c)**a
    w = A1a * dt * r**(-1.0 / a)
    rho, u, P = np.full_like(r, rho0), np.zeros_like(r), np.full_like(r, P0)
    behind = r > R
    near = behind & (w >= tw["w"][0])
    far = behind & ~near
    xi = w[near]**(-a)
    U, C, G = (np.interp(xi, tx["xi"], tx[k]) for k in ("U", "C", "G"))
    rho[near], u[near] = rho0 * G, -r[near] / dt * U
    P[near] = rho[near] * (r[near] / dt)**2 * C / g
    V, Z, G = (np.interp(-w[far], -tw["w"], tw[k]) for k in ("V", "Z", "G"))
    s = A1a * r[far]**(1.0 - 1.0 / a)
    rho[far], u[far], P[far] = rho0 * G, -s * V, rho0 * G * s**2 * Z / g
    return rho, u, P


plotfiles = find_plotfiles(args.plot_dir)
times = np.array([plotfile_time(p) for p in plotfiles])
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.3), constrained_layout=True)
for n, t_want in enumerate(args.times):
    p = plotfiles[int(np.argmin(np.abs(times - t_want)))]
    ds = yt.load(str(p))
    t = float(ds.current_time)
    ray = ds.ortho_ray(0, (1.0e-4, 0.0))
    o = np.argsort(ray[("boxlib", "x")].d)
    x = ray[("boxlib", "x")].d[o]
    sim = [ray[("boxlib", f)].d[o] for f in ("rho", "u", "P")]
    r = np.linspace(1e-3, 1.0, 2000)
    ex = exact(r, t)
    for ax, s_, e_ in zip(axes, sim, ex):
        ax.plot(x, s_, "-", color=colors[n % 5], lw=1.2, label=f"t = {t:.2f}")
        ax.plot(r, e_, "--", color="k", lw=0.8)
for ax, name in zip(axes, (r"density $\rho$", r"velocity $u_r$", r"pressure $P$")):
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("r")
    ax.set_xlim(0, 1)
axes[0].plot([], [], "--", color="k", lw=0.8, label="Guderley")
axes[0].legend(frameon=False, fontsize=8)
fig.savefig(args.out, dpi=200)
print(f"Wrote {args.out}")
