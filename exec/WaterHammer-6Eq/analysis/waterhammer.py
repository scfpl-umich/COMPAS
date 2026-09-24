#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""
Planar water hammer: exact solution and comparison with COMPAS.

Water moving at speed u toward a rigid wall is brought to rest by a planar shock that runs
back into it. For a stiffened gas, p + pinf obeys the ideal-gas Hugoniot, and in the frame
of the incoming water the wall is a piston moving at u. The shock speed relative to the
incoming water is

    W = (gamma+1)/4 u + sqrt( ((gamma+1)/4 u)^2 + c0^2 ),   c0^2 = gamma (p0 + pinf) / rho0,

and mass and momentum across the shock give the water-hammer pressure and density

    p_wh = p0 + rho0 W u,   rho_wh = rho0 W / (W - u).

In the frame of the wall at x = L, the shock leaves the wall at t = 0 and is at
x_s(t) = L - (W - u) t, with the incoming water ahead of it and water at rest at p_wh
behind it.

The script reads the parameters of the run from the configuration log in the plot directory
(or else from the inputs file), the plotfiles (averaged across the strip) and the
waterhammer.csv file that the case writes, prints the exact values next to the simulated
plateau pressure, wall pressure and shock speed, and draws the figure. From the case directory,

    python3 analysis/waterhammer.py plot/WaterHammer-6Eq

writes waterhammer_planar.png (about 4800 px wide) and waterhammer_planar_400dpi.png to the
current directory.
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess as pp                                # noqa: E402  (sets up the Python environment)

import numpy as np                                      # noqa: E402
import matplotlib.pyplot as plt                         # noqa: E402

GRAY = "0.55"
REF = "0.35"      # the gray reference line of the exact solution, as in the Sod figure


# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------

def hugoniot(u, rho0, p0, gamma, pinf):
    """Shock speed W relative to the incoming water, water-hammer pressure and density,
    and the sound speed c0 of the incoming water, for water at (rho0, p0) moving at u
    toward a rigid wall."""
    c0 = np.sqrt(gamma * (p0 + pinf) / rho0)
    a = 0.25 * (gamma + 1.0) * u
    W = a + np.sqrt(a * a + c0 * c0)
    return {"c0": c0, "W": W, "p_wh": p0 + rho0 * W * u, "rho_wh": rho0 * W / (W - u)}


def exact_profile(x, t, L, u, rho0, p0, ex):
    """Pressure, velocity and density along the strip at time t."""
    behind = x > L - (ex["W"] - u) * t
    p = np.where(behind, ex["p_wh"], p0)
    v = np.where(behind, 0.0, u)
    rho = np.where(behind, ex["rho_wh"], rho0)
    return p, v, rho


# ---------------------------------------------------------------------------
# Inputs file and plotfiles
# ---------------------------------------------------------------------------

def read_inputs(path):
    """key = value pairs of an AMReX inputs file, as lists of strings."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                key, val = line.split("=", 1)
                out[key.strip()] = val.replace('"', " ").split()
    return out


def read_run_parameters(plotdir, inputs):
    """Parameters of the run: the Runtime Parameters of the configuration log that COMPAS
    writes to the plot directory, which include any set on the command line (the last value
    of a key wins), or else the inputs file."""
    logs = sorted(glob.glob(os.path.join(plotdir, "config_log_*.txt")))
    if not logs:
        return read_inputs(inputs)
    out, active = {}, False
    with open(logs[-1]) as f:
        for line in f:
            if line.startswith("===="):
                active = "Runtime Parameters" in line
                continue
            m = re.match(r"(\S+)\(nvals = \d+\)\s+::\s+\[(.*)\]", line.strip()) if active else None
            if m:
                out[m.group(1)] = m.group(2).replace(",", " ").split()
    return out


def read_plotfile(path):
    """Time, x and the profiles along x of a 2D plotfile, averaged across the strip."""
    fr = pp.load(path, ("P", "u"))
    return fr.t, fr.x, {k: v.mean(axis=1) for k, v in fr.q.items()}


def read_csv(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}


# ---------------------------------------------------------------------------
# Comparison and figure
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotdir", help="directory with the plotfiles and waterhammer.csv")
    ap.add_argument("--inputs", default="prob/inputs", help="inputs file, used when the plot directory has no configuration log")
    ap.add_argument("--times", type=float, nargs="*", default=None,
                    help="times of the pressure profiles, s (default: five, evenly spaced up to the last plotfile)")
    ap.add_argument("--out", default="waterhammer_planar", help="base name of the figures")
    ap.add_argument("--width-px", type=int, default=4800, help="width of the web figure in pixels")
    args = ap.parse_args()

    inp = read_run_parameters(args.plotdir, args.inputs)
    gamma = float(inp["prob.EOS.gamma_1"][0]); pinf = float(inp["prob.EOS.pinf_1"][0])
    p0 = float(inp["prob.P0"][0]); rho0 = float(inp["prob.rho0"][0]); u = float(inp["prob.u_in"][0])
    L = float(inp["geometry.prob_hi"][0])
    ex = hugoniot(u, rho0, p0, gamma, pinf)
    us = ex["W"] - u

    print("Exact solution (stiffened-gas Hugoniot)")
    print(f"  c0      = {ex['c0']:.2f} m/s")
    print(f"  W       = {ex['W']:.2f} m/s relative to the incoming water, {us:.2f} m/s from the wall")
    print(f"  p_wh    = {ex['p_wh']:.6e} Pa   (acoustic estimate p0 + rho0 c0 u = {p0 + rho0 * ex['c0'] * u:.6e} Pa)")
    print(f"  rho_wh  = {ex['rho_wh']:.2f} kg/m^3")

    # Plotfiles
    snaps = [read_plotfile(p) for _, p in pp.plotfiles(args.plotdir)]
    t_snap = np.array([s[0] for s in snaps])
    times = args.times if args.times else [f * t_snap[-1] for f in (0.2, 0.4, 0.6, 0.8, 1.0)]

    # Plateau behind the shock at the last time, from ten cells behind the exact shock to the wall
    t_end, x, q = snaps[-1]
    dx = x[1] - x[0]
    xs_end = L - us * t_end
    plateau = x > xs_end + 10 * dx
    P_pl = q["P"][plateau]
    print(f"Plateau at t = {t_end * 1e6:.1f} us ({plateau.sum()} cells)")
    print(f"  mean P  = {P_pl.mean():.6e} Pa, relative to p_wh {P_pl.mean() / ex['p_wh'] - 1.0:+.2e}, "
          f"range {P_pl.min() / ex['p_wh'] - 1.0:+.2e} to {P_pl.max() / ex['p_wh'] - 1.0:+.2e}")
    print(f"  mean u  = {q['u'][plateau].mean():+.3e} m/s")

    # Wall pressure and shock trajectory
    csv = read_csv(os.path.join(args.plotdir, "waterhammer.csv"))
    t, xs, Pw = csv["t"], csv["x_shock"], csv["P_wall"]
    late = t >= 0.25 * t[-1]
    print(f"Wall pressure for t >= {0.25 * t[-1] * 1e6:.0f} us")
    print(f"  mean    = {Pw[late].mean():.6e} Pa, relative to p_wh {Pw[late].mean() / ex['p_wh'] - 1.0:+.2e}")
    slope, x0 = np.polyfit(t[late], xs[late], 1)
    print(f"Shock trajectory for t >= {0.25 * t[-1] * 1e6:.0f} us")
    print(f"  speed from the wall = {-slope:.2f} m/s, relative to the exact {-slope / us - 1.0:+.2e}")
    print(f"  W = {-slope + u:.2f} m/s, relative to the exact {(-slope + u) / ex['W'] - 1.0:+.2e}")
    print(f"  position at t = {t[-1] * 1e6:.1f} us: {xs[-1]:.4f} m, exact {L - us * t[-1]:.4f} m "
          f"({(xs[-1] - (L - us * t[-1])) / dx:+.2f} cells)")

    # Figure: pressure along the strip, each time labeled at its shock, and the wall pressure
    pp.use_style()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.28})
    GPa = 1e-9
    p_ac = p0 + rho0 * ex["c0"] * u   # acoustic estimate, for reference
    top = max(1.27 * ex["p_wh"], 1.05 * Pw.max()) * GPa

    ax1.axhline(ex["p_wh"] * GPa, color=GRAY, lw=0.8, ls=(0, (4, 3)), zorder=0)
    ax1.text(0.012 * L, ex["p_wh"] * GPa + 0.02 * top, r"$p_{wh}$", color="0.35", va="bottom")
    xe = np.linspace(0.0, L, 4001)
    for n, tt in enumerate(times):
        j = int(np.argmin(np.abs(t_snap - tt)))
        ts, xx, qq = snaps[j]
        c = colors[n % len(colors)]
        pe, _, _ = exact_profile(xe, ts, L, u, rho0, p0, ex)
        ax1.plot(xe, pe * GPa, color=REF, lw=1.2, ls=(0, (3, 2)), zorder=3)
        ax1.plot(xx, qq["P"] * GPa, color=c, lw=1.4, zorder=2)
        ax1.text(L - us * ts - 0.013 * L, 0.49 * ex["p_wh"] * GPa, rf"${ts * 1e6:.0f}\,\mu$s", color=c, ha="right", va="center")
    ax1.plot([], [], color="0.35", lw=1.4, label="COMPAS")
    ax1.plot([], [], color=REF, lw=1.2, ls=(0, (3, 2)), label="exact")
    ax1.set_xlim(0.0, L)
    ax1.set_ylim(-0.035 * top, top)
    ax1.set_xlabel(r"$x$ (m)")
    ax1.set_ylabel(r"$p$ (GPa)")
    ax1.legend(loc="upper right", ncol=2, columnspacing=1.2, handlelength=1.8)
    ax1.set_title("(a)", loc="left")

    ax2.axhline(p_ac * GPa, color=GRAY, lw=0.8, ls=(0, (1, 2)), zorder=0)
    ax2.plot(t * 1e6, Pw * GPa, color=colors[0], lw=1.4, zorder=2, label="COMPAS")
    ax2.axhline(ex["p_wh"] * GPa, color=REF, lw=1.2, ls=(0, (3, 2)), zorder=3, label=r"$p_{wh}$, exact")
    ax2.plot([], [], color=GRAY, lw=0.8, ls=(0, (1, 2)), label=r"$p_0 + \rho_0 c_0 u$, acoustic")
    ax2.set_xlim(0.0, t[-1] * 1e6)
    ax2.set_ylim(-0.035 * top, top)
    ax2.set_xlabel(r"$t$ ($\mu$s)")
    ax2.set_ylabel(r"$p_{\rm wall}$ (GPa)")
    ax2.legend(loc="lower right")
    ax2.set_title("(b)", loc="left")

    pp.save(fig, f"{args.out}.png", args.width_px)


if __name__ == "__main__":
    main()
