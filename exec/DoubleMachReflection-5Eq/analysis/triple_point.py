#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Primary triple point and foot of the Mach stem at t = 0.2, against three-shock theory, and the
close-up of the jet, the slip lines and the Mach stem.

In each row of the finest grid, the leading shock is the furthest point in x where the
pressure crosses halfway between the pressures ahead of and behind the incident shock. Above
the primary triple point this is the incident shock, whose exact position is
x_s(y, t) = x0 + y tan(theta_w) + Us t / cos(theta_w); below it, the Mach stem runs ahead.
The offset d(y) of the leading shock from x_s is therefore a constant (a fraction of a cell)
above the triple point and grows below it. The triple point is where a straight line fitted
to d(y) just below the kink meets the constant, and the foot of the Mach stem is the leading
shock in the first row. The trajectory angle of the triple point, seen from the wedge tip, is
chi = atan(y_T / (x_T - x0)). The spread is over kinks placed 2, 3 and 4 cells ahead.

Three-shock theory (von Neumann; see Ben-Dor, Shock Wave Reflection Phenomena, 2007) gives chi
when the Mach stem is taken as straight and normal to the wall: in the frame of the triple
point, the gas ahead has the Mach number M0 = Ms / cos(theta_w + chi), meets the incident shock
at phi_1 = 90 - theta_w - chi degrees and the Mach stem at phi_3 = 90 - chi degrees, and chi is
the angle for which the pressures and the flow directions behind the reflected shock and the
Mach stem agree. The Mach stem is curved in the simulations, so this is a model, not an exact
solution.

With --figure, one panel per run over [2.0, 2.9] x [0, 0.5], in rows of --ncols panels (default
all in one row): a light gray schlieren of the density (|grad rho|, logarithmic, the same limits in
every panel) under 30 density contours from 1.73 to 21.5, the measured triple point and foot of the
Mach stem (open circles) and the undisturbed incident shock x_s(y, t) below the triple point
(dashed). Each panel is labeled with the finest cell size of its run (--labels), and --row-labels
names each row, for example by its scheme.

    python3 analysis/triple_point.py plot/DoubleMach240 plot/DoubleMach \
        plot/DoubleMachWENO5_240 plot/DoubleMachWENO5 --figure dmr_zoom.png --ncols 2 \
        --labels 1/240 1/480 1/240 1/480 --row-labels "MUSCL-MC, HLLC" "WENO5, HLL"
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess as pp                                # noqa: E402  (sets up the Python environment)

import numpy as np                                      # noqa: E402
from scipy.optimize import brentq                       # noqa: E402

GAMMA, MS, THETA_W, X0, P0, RHO0 = 1.4, 10.0, 30.0, 1.0 / 6.0, 1.0, 1.4
RHO_LEVELS = np.linspace(1.73, 21.5, 30)                # as in Woodward and Colella (1984)


def oblique(M, phi, g=GAMMA):
    """Pressure ratio, flow deflection and downstream Mach number of an oblique shock with
    upstream Mach number M at the angle phi (radians) to the flow."""
    Mn = M * np.sin(phi)
    pr = 1.0 + 2.0 * g / (g + 1.0) * (Mn**2 - 1.0)
    th = np.arctan(2.0 / np.tan(phi) * (Mn**2 - 1.0) / (M**2 * (g + np.cos(2.0 * phi)) + 2.0))
    Mn2 = np.sqrt((1.0 + 0.5 * (g - 1.0) * Mn**2) / (g * Mn**2 - 0.5 * (g - 1.0)))
    return pr, th, Mn2 / np.sin(phi - th)


def three_shock_chi(Ms=MS, theta_w=THETA_W):
    """Triple-point trajectory angle (degrees) from three-shock theory with a straight Mach
    stem normal to the wall."""
    tw = np.radians(theta_w)

    def mismatch(chi):
        M0 = Ms / np.cos(tw + chi)
        p1, th1, M1 = oblique(M0, 0.5 * np.pi - tw - chi)         # incident shock
        p3, th3, _ = oblique(M0, 0.5 * np.pi - chi)                # Mach stem
        # reflected shock: the weakest one that brings the pressure up to p3
        f = lambda phi: p1 * oblique(M1, phi)[0] - p3
        phi = np.linspace(np.arcsin(1.0 / M1) + 1e-9, 0.5 * np.pi, 4001)
        k = np.where(np.diff(np.sign(f(phi))) != 0)[0][0]
        th2 = oblique(M1, brentq(f, phi[k], phi[k + 1]))[1]
        return th1 - th2 - th3                                     # flow directions agree

    return np.degrees(brentq(mismatch, np.radians(1.0), np.radians(30.0)))


def incident_shock(y, t):
    """x of the undisturbed incident shock."""
    tw = np.radians(THETA_W)
    return X0 + y * np.tan(tw) + MS * np.sqrt(GAMMA * P0 / RHO0) * t / np.cos(tw)


def leading_shock(x, P, P_thr):
    """Furthest x in each row where P crosses P_thr, by linear interpolation (P is [i, j])."""
    out = np.full(P.shape[1], np.nan)
    for j in range(P.shape[1]):
        i = np.where(P[:, j] > P_thr)[0]
        if len(i) and i[-1] + 1 < len(x):
            i = i[-1]
            out[j] = x[i] + (x[i + 1] - x[i]) * (P[i, j] - P_thr) / (P[i, j] - P[i + 1, j])
    return out


def measure(fr, n_kink=3.0):
    """(x_foot, x_T, y_T, chi) in a frame with the pressure. The triple point is first placed
    where d(y) exceeds its value on the undisturbed shock by n_kink cells."""
    x, y, dx = fr.x, fr.y, fr.dx
    PS = P0 * (2.0 * GAMMA * MS**2 - (GAMMA - 1.0)) / (GAMMA + 1.0)
    xs = incident_shock(y, fr.t)
    off = leading_shock(x, fr.q["P"], 0.5 * (P0 + PS)) - xs
    d_far = np.nanmedian(off[(y > 0.65) & (y < 0.95)])        # well above the triple point
    y0 = y[np.where(off > d_far + n_kink * dx)[0].max()]       # highest row where the stem leads
    sel = (y > y0 - 0.08) & (y < y0 - 0.01) & np.isfinite(off)
    a, b = np.polyfit(y[sel], off[sel], 1)                     # d(y) just below the kink
    y_T = (d_far - b) / a
    x_T = incident_shock(y_T, fr.t) + d_far
    x_foot = xs[0] + off[0]
    return x_foot, x_T, y_T, np.degrees(np.arctan2(y_T, x_T - X0))


def zoom_figure(runs, labels, schemes, window, out, ncols=None, row_labels=None):
    x0w, x1w, y0w, y1w = window
    nc = ncols or len(runs)
    nr = -(-len(runs) // nc)
    pp.use_style()
    fig, axes = pp.panel_grid(nr, nc, (y1w - y0w) / (x1w - x0w), left=0.55 + (0.3 if row_labels else 0.0),
                              right=0.1, top=0.1, bottom=0.5, wgap=0.12, hgap=0.12)
    for k, ax in enumerate(axes.flat):
        r, c = divmod(k, nc)
        if k >= len(runs):
            ax.set_visible(False)
            continue
        fr, m = runs[k]
        ix = (fr.x > x0w - 0.02) & (fr.x < x1w + 0.02)
        jy = fr.y < y1w + 0.02
        rho, xs, ys = fr.q["rho"][np.ix_(ix, jy)], fr.x[ix], fr.y[jy]
        shade = pp.schlieren(rho, fr.dx, lims=(8.0, 3000.0), darkest=0.5)
        ax.imshow(1.0 - shade.T, cmap="gray", vmin=0.0, vmax=1.0, origin="lower", interpolation="antialiased",
                  interpolation_stage="data", extent=[xs[0] - fr.dx / 2, xs[-1] + fr.dx / 2, ys[0] - fr.dx / 2,
                                                       ys[-1] + fr.dx / 2])
        ax.contour(xs, ys, rho.T, levels=RHO_LEVELS, colors=pp.INK, linewidths=0.35)
        x_foot, x_T, y_T, _ = m
        yy = np.array([0.0, y_T])
        ax.plot(incident_shock(yy, fr.t), yy, color="0.55", lw=0.7, ls=(0, (4, 3)), zorder=4)
        ax.plot([x_T, x_foot], [y_T, 0.004], "o", ms=5, mfc="none", mec="#E1701A", mew=0.9, zorder=5, clip_on=False)
        label = rf"$\Delta x = {labels[k] if labels else f'{fr.dx:.5f}'}$"
        if schemes:
            label = f"{schemes[k]},  " + label
        ax.text(0.03, 0.97, label, transform=ax.transAxes, ha="left", va="top", fontsize=10, color=pp.INK,
                bbox=dict(fc="white", ec="none", pad=1.5, alpha=0.85), zorder=20)
        ax.set_xlim(x0w, x1w)
        ax.set_ylim(y0w, y1w)
        ax.set_xticks(np.arange(x0w, x1w - 0.05, 0.2))           # no label at the shared edge
        if r > 0:
            ax.set_yticks(np.arange(y0w, y1w - 0.05, 0.1))       # nor at the edge shared with the row above
        ax.set_aspect("equal")
        if r == nr - 1 or k + nc >= len(runs):
            ax.set_xlabel(r"$x$")
        else:
            ax.tick_params(labelbottom=False)
        if c == 0:
            ax.set_ylabel(r"$y$")
        else:
            ax.tick_params(labelleft=False)
    for r, text in enumerate((row_labels or [])[:nr]):
        pos = axes[r, 0].get_position()
        fig.text(0.1 / fig.get_figwidth(), 0.5 * (pos.y0 + pos.y1), text, rotation=90, ha="left", va="center",
                 fontsize=10, color=pp.INK)
    pp.save(fig, out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="plot directories, e.g. plot/DoubleMach")
    ap.add_argument("--figure", help="write the close-up figure to this file")
    ap.add_argument("--labels", nargs="+", help="finest dx of each run, e.g. 1/240 1/480")
    ap.add_argument("--schemes", nargs="+", help="scheme of each run, written in its panel, e.g. 'WENO5, HLL'")
    ap.add_argument("--ncols", type=int, help="panels per row of the figure (default all in one row)")
    ap.add_argument("--row-labels", nargs="+", help="label to the left of each row of panels, e.g. its scheme")
    ap.add_argument("--window", type=float, nargs=4, default=[2.0, 2.9, 0.0, 0.5], help="x0 x1 y0 y1")
    args = ap.parse_args()

    chi_3st = three_shock_chi()
    tan_c = np.tan(np.radians(chi_3st))
    L = incident_shock(0.0, 0.2) - X0
    y_3st = L * tan_c / (1.0 - np.tan(np.radians(THETA_W)) * tan_c)
    print(f"three-shock theory (straight Mach stem normal to the wall): chi = {chi_3st:.2f} deg, "
          f"which puts the triple point at x = {X0 + y_3st / tan_c:.3f}, y = {y_3st:.3f} at t = 0.2")
    runs = []
    for p in args.runs:
        fr = pp.load(pp.pick(p, 0.2)[1], ("rho", "P"))
        rows = np.array([measure(fr, n) for n in (2.0, 3.0, 4.0)])
        x_foot, x_T, y_T, chi = rows[1]
        spread = np.ptp(rows[:, 1:], axis=0)
        print(f"{p}: t = {fr.t:.4f}, dx = 1/{1 / fr.dx:.0f}: Mach stem foot x = {x_foot:.4f}, "
              f"triple point x = {x_T:.4f}, y = {y_T:.4f} (spread {spread[0]:.4f}, {spread[1]:.4f}), "
              f"chi = {chi:.2f} deg (spread {spread[2]:.2f})")
        runs.append((fr, rows[1]))
    if args.figure:
        zoom_figure(runs, args.labels, args.schemes, args.window, args.figure, args.ncols, args.row_labels)
