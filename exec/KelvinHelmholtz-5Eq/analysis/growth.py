#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Growth rate of the seeded mode from the in-situ history mode.csv of each run, and the figure
of the growth rate against Blumen (1970).

mode.csv has one line every run.user_output_int coarse steps, with the amplitude A of the seeded
Fourier mode of v (integrated over y) and the largest rho v^2/2 (see prob/ProblemICBC.H). The
growth rate is the least-squares slope of ln A over the linear phase, taken from the first time
A reaches 100 A(0), after the acoustic transient of the seed (which is not the eigenmode) has
died out, to the first time the largest |v| reaches 0.01 U, before the mode becomes nonlinear.
The uncertainty is the largest of three estimates: the standard error of the slope; half the
difference between the slopes of the first and second halves of the window, which measures how
far ln A is from a straight line; and half the spread of the slopes when the window starts at
30, 100 or 300 A(0) and ends at max |v| = 0.005, 0.01 or 0.02 U, which measures how much the
residual ringing of the seed and the onset of nonlinearity bias the fit.

Each run is given as label:M:dx:path, with dx given as delta/dx at the finest level. The rates
are written to --out. With --figure, the figure has
(a) the growth rate against M: Blumen's Table 1, the maximum growth rate of his pressure
    equation (7) solved by blumen.py, and the rates with their uncertainty; below, the relative
    difference from the rate of equation (7) at the simulated wavenumber, in percent, with the
    band of the three-digit rounding of Table 1;
(b) the amplitude of the seeded mode against time for the run labeled --example, with the fit
    window shaded and the fitted exponential.
The panels carry only the letters (a) and (b); the caption says what each shows.

    python3 analysis/growth.py 0.2:0.2:16.27:plot/KH-M0.2/mode.csv 0.4:0.4:16.14:plot/KH-M0.4/mode.csv ... \
        --out growth_rates.csv --figure kh_growth.png --example 0.4
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess as pp                                # noqa: E402  (sets up the Python environment)

import numpy as np                                      # noqa: E402
import matplotlib.pyplot as plt                         # noqa: E402

import blumen                                           # noqa: E402

A_START = 100.0   # the window starts where A = A_START A(0)
V_END = 0.01      # and ends where max |v| = V_END U
A_STARTS = (30.0, 100.0, 300.0)    # other starts and ends, for the sensitivity of the rate to the window
V_ENDS = (0.005, 0.01, 0.02)


def slope(t, y):
    """Least-squares slope of y(t) and its standard error."""
    p, cov = np.polyfit(t, y, 1, cov=True)
    return p[0], np.sqrt(cov[0, 0])


def window(t, A, vmax, a_start, v_end):
    """Indices [i0, i1) of the fit window, or None if the run does not cover it."""
    if not (A >= a_start * A[0]).any():
        return None
    i0 = np.argmax(A >= a_start * A[0])
    i1 = np.argmax(vmax >= v_end) if (vmax >= v_end).any() else len(t)
    return (i0, i1) if i1 - i0 >= 10 else None


def fit(csv):
    """Growth rate, its uncertainty and the fit window (t0, t1) for one mode.csv."""
    d = np.genfromtxt(csv, delimiter=",", names=True)
    t, A = d["t"], d["A"]
    vmax = np.sqrt(2.0 * d["Ey_max"])          # rho = 1 to O(M^2 eps^2)
    w = window(t, A, vmax, A_START, V_END)
    if w is None:
        raise SystemExit(f"{csv}: the run does not cover the fit window")
    tw, lw = t[w[0]:w[1]], np.log(A[w[0]:w[1]])
    s, se = slope(tw, lw)
    h = len(tw) // 2
    s1, _ = slope(tw[:h], lw[:h])
    s2, _ = slope(tw[h:], lw[h:])
    others = [window(t, A, vmax, a, v) for a in A_STARTS for v in V_ENDS]
    others = [slope(t[i0:i1], np.log(A[i0:i1]))[0] for i0, i1 in filter(None, others)]
    return s, max(se, 0.5 * abs(s2 - s1), 0.5 * np.ptp(others)), tw[0], tw[-1]


def figure(r, mode_csv, k, out):
    """The growth-rate figure of the rates r, with the run of row k and its mode.csv in panel (b)."""
    blue, orange, ref = "#4477AA", "#E1701A", "0.55"
    table = blumen.TABLE1

    # Maximum growth rate of equation (7), and its value at the wavenumber of each run (Blumen's alpha_max)
    M_curve = np.linspace(0.0, 0.95, 20)
    s_curve = np.array([blumen.most_unstable(m)[1] for m in M_curve])
    alpha_run = np.interp(r["M"], table[:, 0], table[:, 1])
    s_theory = np.array([blumen.growth_rate(a, m) for a, m in zip(alpha_run, r["M"])])
    s_table = np.interp(r["M"], table[:, 0], table[:, 2])
    for row, st, sb in zip(r, s_theory, s_table):
        print(f"M = {row['M']:.1f}, delta/dx = {row['delta_over_dx']:.2f}: COMPAS {row['rate']:.5f} +- "
              f"{row['rate_err']:.5f}, equation (7) {st:.5f}, Table 1 {sb:.3f}, "
              f"difference {100 * (row['rate'] / st - 1):+.2f} %")

    pp.use_style()
    fig = plt.figure(figsize=(10.0, 4.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.15], width_ratios=[1, 1], hspace=0.08, wspace=0.26)
    ax1 = fig.add_subplot(gs[0, 0])
    ax3 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax2 = fig.add_subplot(gs[:, 1])

    # (a) growth rate against M
    coarse = r["delta_over_dx"] < 24
    fine = ~coarse
    ax1.plot(M_curve, s_curve, color=ref, lw=1.0, label="Blumen's equation (7), solved here")
    ax1.plot(table[:, 0], table[:, 2], "o", mfc="none", mec="k", ms=7, mew=0.9, zorder=4,
             label="Blumen (1970), Table 1")
    ax1.errorbar(r["M"][coarse], r["rate"][coarse], yerr=r["rate_err"][coarse], fmt="s", color=blue, ms=3.8,
                 label=r"COMPAS, $\Delta x = \delta/16$")
    ax1.errorbar(r["M"][fine] + 0.025, r["rate"][fine], yerr=r["rate_err"][fine], fmt="D", color=orange, ms=3.8,
                 label=r"COMPAS, $\Delta x = \delta/32$")
    ax1.set_xlim(-0.03, 1.03)
    ax1.set_ylim(-0.009, 0.21)                     # room for the Table 1 point at M = 1
    ax1.set_ylabel(r"growth rate, $\alpha c_i\ (U/\delta)$")
    ax1.legend(loc="lower left")
    ax1.tick_params(labelbottom=False)
    ax1.set_title("(a)", loc="left")

    # Difference from equation (7) at the simulated wavenumber, with the rounding of Table 1
    M_band = np.linspace(0.1, 0.9, 100)
    half = 100 * 0.0005 / np.interp(M_band, table[:, 0], table[:, 2])
    ax3.fill_between(M_band, -half, half, color="0.9", lw=0, label="rounding of Table 1")
    ax3.axhline(0.0, color=ref, lw=0.8)
    d = 100 * (r["rate"] / s_theory - 1)
    e = 100 * r["rate_err"] / s_theory
    ax3.errorbar(r["M"][coarse], d[coarse], yerr=e[coarse], fmt="s", color=blue, ms=4)
    ax3.errorbar(r["M"][fine] + 0.025, d[fine], yerr=e[fine], fmt="D", color=orange, ms=3.5)
    ax3.set_ylim(-1.6, 1.6)
    ax3.set_xlabel(r"Mach number, $M = U/a$")
    ax3.set_ylabel("difference\nfrom (7), %")
    ax3.legend(loc="lower left", fontsize=8, handlelength=1.2)

    # (b) amplitude of the seeded mode against time
    m = np.genfromtxt(mode_csv, delimiter=",", names=True)
    M, rate, t0, t1 = float(r["M"][k]), float(r["rate"][k]), float(r["t0"][k]), float(r["t1"][k])
    rate_7 = float(s_theory[k])
    res = 16 if coarse[k] else 32
    t, A = m["t"], m["A"] / m["A"][0]
    w = (t >= t0) & (t <= t1)
    c = np.exp(np.mean(np.log(A[w]) - rate * t[w]))
    ax2.axvspan(t0, t1, color="0.92", lw=0)          # the fit window
    ax2.semilogy(t, A, color=blue, lw=1.6, label=rf"COMPAS, $M = {M:g}$, $\Delta x = \delta/{res}$")
    tt = np.linspace(0, t[-1], 200)
    ax2.semilogy(tt, c * np.exp(rate * tt), color="k", lw=1.0, ls=(0, (7, 3)),
                 label=rf"fit over the window, $\alpha c_i$ = {rate:.4f}" + "\n" + rf"(equation (7): {rate_7:.4f})")
    ax2.set_xlim(0, t[-1])
    ax2.set_ylim(0.5, 3e5)
    ax2.set_xlabel(r"time, $t\ (\delta/U)$")
    ax2.set_ylabel(r"mode amplitude, $A(t)/A(0)$")
    ax2.legend(loc="upper left")
    ax2.set_title("(b)", loc="left")
    pp.save(fig, out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label:M:dx:path/to/mode.csv, dx = delta/dx at the finest level")
    ap.add_argument("--out", default="growth_rates.csv", help="growth rates (CSV)")
    ap.add_argument("--figure", help="write the growth-rate figure to this file")
    ap.add_argument("--example", default="0.4", help="label of the run shown in panel (b)")
    args = ap.parse_args()

    rows, paths = [], []
    for run in args.runs:
        label, M, dx, path = run.split(":", 3)
        s, e, t0, t1 = fit(path)
        rows.append((label, float(M), float(dx), s, e, t0, t1))
        paths.append(path)
        print(f"{label:>16s}  M = {float(M):.1f}  dx = delta/{dx}  rate {s:.5f} +- {e:.5f}  (t = {t0:.1f} to {t1:.1f})")
    with open(args.out, "w") as f:
        f.write("label,M,delta_over_dx,rate,rate_err,t0,t1\n")
        for row in rows:
            f.write("{},{},{},{:.6f},{:.6f},{:.3f},{:.3f}\n".format(*row))
    print("wrote", args.out)
    if args.figure:
        labels = [row[0] for row in rows]
        if args.example not in labels:
            raise SystemExit(f"--example {args.example}: no run with this label")
        k = labels.index(args.example)
        r = np.atleast_1d(np.genfromtxt(args.out, delimiter=",", names=True, dtype=None, encoding=None))
        figure(r, paths[k], k, args.figure)
