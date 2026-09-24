#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Measurements and figures of the shock refraction runs, one run per interface angle beta in
<root>/b<beta>/plot, with the plotfiles and the waves.csv that the case writes every coarse step
(run.user_output_int = 1). The interface runs from the apex of the wedge at the origin at the
angle 90 deg - beta to the wall; s is the distance along it from the apex and n the distance from
it into the air.

waves <root>     measures the wave pattern of each run, prints it, writes it to
                 <root>/measurements.json (or --json) and draws its figure (--out):
  * Incident shock speed, from x_shock in waves.csv (the top row of cells), against U_s.
  * Precursor: the slope of S_t(t) in waves.csv, the front of the transmitted wave along the bottom
    wall, fitted while it is between 10 and 70 mm (the domain ends at 80 mm), against the sound
    speed of the water.
  * Wave angles at the last plotfile, on arcs of radius r around the refraction point R, where the
    incident shock (the top row) meets the interface. On each arc the reflected shock in the air
    (and in RRR the transmitted shock in the water) is where the pressure crosses halfway between
    the states on either side of it. A straight line fitted to these crossings gives the angle of
    the wave to the interface: omega_r for the reflected shock and phi_b for the transmitted one.
  * Front of the transmitted wave 2 mm under the interface, where the pressure first exceeds
    P0 + 0.05 (P1 - P0), in each plotfile: its speed along the interface is the speed of the
    refraction point, U_s / sin(beta), when the wave is bound to R (RRR, BPR), and the water sound
    speed when it is a free precursor.
  * Wall pressure levels: in each plotfile, the end of the stretch of the bottom wall, continuous
    from the apex, where the pressure exceeds P0 + f (P1 - P0), for f = 0.01, 0.05, 0.25 and 0.5,
    and the speed of each level fitted while it is between 10 and 70 mm. The 1 percent level is
    the leading edge of the precursor.
  * Triple point: in columns a distance d behind the incident shock, the highest point where the
    pressure rises above 1.25 P1 (P1 the pressure behind the incident shock) is on the reflected
    shock; a straight line through these points, extrapolated to d = 0, meets the incident shock at
    the triple point T. The trajectory angle chi is the angle of T above the interface seen from
    the apex, and the slope of a line fitted to T(t) in interface coordinates over the second half
    of the run gives it without assuming the apex as origin. In the regular reflections the same
    construction gives a point about 0.4 mm (three cells) off the interface, the offset of the
    method.

  The figure:
  (a) wave angles to the interface against beta, with polar theory (polar.py) and the RRR to BPR
      and detachment (FPR to FMR) transitions;
  (b) the leading edge of the precursor along the bottom wall, with a wave leaving the apex at the
      water sound speed when the incident shock reaches it (t0);
  (c) the triple point of the beta = 50 deg run, with a line fitted over the second half of the
      run and chi = 2 +/- 0.1 deg, the trajectory angle measured in the experiment, drawn over the
      fitted stretch from the start of the fitted line.

snapshots <root> the four runs at their last plotfile: the air volume fraction through the lighter
  teal-to-sand map with its 0.5 contour, darkened by a light schlieren of |grad rho| in the pure air
  and of |grad P| elsewhere, in the water and across the interface, where a density schlieren would
  saturate, with fixed limits in every panel, and every AMR patch box, in the color of its level,
  on the beta = 50 deg panel. Each panel is labeled with beta, the pattern and t - t0.

    python3 analysis/refraction.py waves plot --out shockrefraction_waves.png
    python3 analysis/refraction.py snapshots plot --out shockrefraction_snapshots.png
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess as pp                                # noqa: E402  (sets up the Python environment)

import numpy as np                                      # noqa: E402
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.colors import Normalize                 # noqa: E402
from scipy.interpolate import RegularGridInterpolator   # noqa: E402
from scipy.ndimage import minimum_filter                # noqa: E402

import polar                                            # noqa: E402

P0 = polar.P0
U_S = polar.MS * polar.sound_speed(polar.AIR)                           # incident shock speed, m/s
P1 = P0 * (1.0 + 2.8 / 2.4 * (polar.MS**2 - 1.0))                       # pressure behind it
C_WATER = polar.sound_speed(polar.WATER)                                # water sound speed, m/s
T0 = 2.0e-3 / U_S           # the shock starts 2 mm ahead of the apex (prob.x_s) and reaches it at T0
WALL_LEVELS = (0.01, 0.05, 0.25, 0.5)   # fractions of P1 - P0 followed along the bottom wall
RUNS = [(15.0, "RRR"), (20.5, "FPR"), (30.0, "FPR"), (50.0, "FMR")]


def run_dir(root, beta):
    return Path(root) / f"b{beta:g}" / "plot"


# ---------------------------------------------------------------------------
# Measurements
# ---------------------------------------------------------------------------

def to_interface(x, y, beta):
    """Distance s along the interface from the apex, and n from the interface into the air."""
    tw = np.radians(90.0 - beta)
    return x * np.cos(tw) + y * np.sin(tw), -x * np.sin(tw) + y * np.cos(tw)


def from_interface(s, n, beta):
    tw = np.radians(90.0 - beta)
    return s * np.cos(tw) - n * np.sin(tw), s * np.sin(tw) + n * np.cos(tw)


def crossing(x, p, level):
    """Largest x where p is above level, by linear interpolation to the next point (nan if none)."""
    i = np.where(p > level)[0]
    if len(i) == 0 or i[-1] + 1 >= len(x):
        return np.nan
    i = i[-1]
    return x[i] + (x[i + 1] - x[i]) * (p[i] - level) / (p[i] - p[i + 1])


def wall_front(x, p, level):
    """End of the stretch of the wall, continuous from the apex (x = 0), where p is above level, by
    linear interpolation to the next point (nan if there is none, or if it reaches the end)."""
    i0 = int(np.searchsorted(x, 0.0))
    above = p[i0:] > level
    if not above[0] or above.all():
        return np.nan
    i = i0 + int(np.argmin(above)) - 1
    return x[i] + (x[i + 1] - x[i]) * (p[i] - level) / (p[i] - p[i + 1])


def fit_line(t, x):
    """Slope, its standard error and the intercept of a least-squares line x = a + b t."""
    A = np.vstack([np.ones_like(t), t]).T
    c, *_ = np.linalg.lstsq(A, x, rcond=None)
    r = x - A @ c
    cov = np.linalg.inv(A.T @ A) * (r @ r) / max(len(t) - 2, 1)
    return c[1], np.sqrt(cov[1, 1]), c[0]


def speeds(csv):
    """Incident shock speed, and precursor speed along the bottom wall with the spread of the
    speeds fitted to the first and second halves of the window."""
    d = np.genfromtxt(csv, delimiter=",", names=True)
    out = {}
    m = (d["x_shock"] > -1e29) & (d["t"] > 2e-6)
    out["U_s"], out["U_s_err"], _ = fit_line(d["t"][m], d["x_shock"][m])
    S = d["S_t"]
    m = (S > 0.010) & (S < 0.070)
    if m.sum() > 10:
        t, s = d["t"][m], S[m]
        v, e, a = fit_line(t, s)
        h = len(t) // 2
        v1, v2 = fit_line(t[:h], s[:h])[0], fit_line(t[h:], s[h:])[0]
        out.update(S_t_speed=v, S_t_err=max(e, 0.5 * abs(v1 - v2)), S_t_speed_halves=[v1, v2],
                   S_t_origin=-a / v)   # time at which the fitted front leaves the apex, x = 0
    return out


def arc_crossings(interp, xR, yR, beta, side, radii):
    """Points (s - s_R, n) where the pressure on arcs around R crosses halfway across the
    reflected shock (side = +1, air) or the transmitted shock (side = -1, water)."""
    tw = np.radians(90.0 - beta)
    t_hat = np.array([np.cos(tw), np.sin(tw)])      # along the interface, away from the apex
    n_hat = side * np.array([-np.sin(tw), np.cos(tw)])
    psi_max = (180.0 - beta - 3.0) if side > 0 else 178.0
    psi = np.radians(np.linspace(1.0, psi_max, 4000))
    pts = []
    for r in radii:
        xy = np.array([xR, yR])[:, None] + r * (-np.outer(t_hat, np.cos(psi)) + np.outer(n_hat, np.sin(psi)))
        p = interp(xy.T)
        k = int(np.argmax(p))
        # from behind the reflected shock to the incident-shock state, or behind the transmitted
        # shock to the water at rest
        lo = np.median(p[psi > psi[-1] - np.radians(10.0)]) if side > 0 else P0
        mid = 0.5 * (p[k] + lo)
        j = np.where(p[k:] < mid)[0]
        if len(j) == 0:
            continue
        j = k + j[0]
        w = (p[j - 1] - mid) / (p[j - 1] - p[j])
        ps = psi[j - 1] + w * (psi[j] - psi[j - 1])
        pts.append((-r * np.cos(ps), r * np.sin(ps)))
    return np.array(pts)


def wave_angle(pts):
    """Angle (deg) to the interface of the line fitted to the arc crossings, and its uncertainty
    (the larger of the standard error and half the difference between the inner and outer halves)."""
    b, e, _ = fit_line(-pts[:, 0], pts[:, 1])
    h = len(pts) // 2
    b1, b2 = fit_line(-pts[:h, 0], pts[:h, 1])[0], fit_line(-pts[h:, 0], pts[h:, 1])[0]
    err = max(np.degrees(e / (1 + b * b)), 0.5 * abs(np.degrees(np.arctan(b1) - np.arctan(b2))))
    return np.degrees(np.arctan(b)), err


def triple_point(x, y, P, beta, xi, levels=(1.25,)):
    """Triple point (x_T, y_T) for each threshold level (in units of P1), from the reflected shock
    in columns 0.5 to 2 mm behind the incident shock at xi."""
    dx = x[1] - x[0]
    out = []
    for lev in levels:
        ds, ys = [], []
        for d in np.arange(0.5e-3, 2.0e-3 + 1e-9, dx):
            i = int(np.argmin(np.abs(x - (xi - d))))
            col = P[i]
            m = np.where((col > lev * P1) & (y > x[i] * np.tan(np.radians(90.0 - beta))))[0]  # above the interface
            if len(m) == 0 or m[-1] + 1 >= len(y):
                continue
            j = m[-1]
            ds.append(xi - x[i])
            ys.append(y[j] + (y[j + 1] - y[j]) * (col[j] - lev * P1) / (col[j] - col[j + 1]))
        out.append((xi, fit_line(np.array(ds), np.array(ys))[2]))
    return out


def measure(root, betas, out):
    res = {"U_s_theory": U_S, "c_water": C_WATER, "T0": T0, "runs": {}}
    for beta in betas:
        tag = f"{beta:g}"
        run = run_dir(root, beta)
        r = speeds(run / "waves.csv")
        pf = pp.plotfiles(run)

        # Wave angles at the last plotfile
        fr = pp.load(pf[-1][1], ("P",))
        x, y, P = fr.x, fr.y, fr.q["P"]
        interp = RegularGridInterpolator((x, y), P, bounds_error=False, fill_value=np.nan)
        xi = crossing(x, P[:, -1], 0.5 * (P0 + P1))          # incident shock along the top row
        xR, yR = xi, xi * np.tan(np.radians(90.0 - beta))
        r["t_last"] = fr.t
        radii = np.linspace(1.0e-3, 8.0e-3, 29)
        if beta < 41.0:     # regular reflection in the air
            r["omega_r"], r["omega_r_err"] = wave_angle(arc_crossings(interp, xR, yR, beta, +1, radii))
        if beta < 19.7:     # transmitted shock
            r["phi_b"], r["phi_b_err"] = wave_angle(arc_crossings(interp, xR, yR, beta, -1, radii))

        # Transmitted front and wall levels in every plotfile once the pattern is 5 us old, and the
        # triple point once it is 10 us old, with the spread over the thresholds
        traj, fronts = [], []
        for tp, f in pf:
            if tp < T0 + 5e-6:
                continue
            fr = pp.load(f, ("P",))
            x, y, P = fr.x, fr.y, fr.q["P"]
            interp = RegularGridInterpolator((x, y), P, bounds_error=False, fill_value=np.nan)
            s_line = np.linspace(0.0, 0.12, 12001)
            p_line = interp(np.c_[from_interface(s_line, -2e-3 + 0 * s_line, beta)])
            fronts.append([tp, crossing(s_line, p_line, P0 + 0.05 * (P1 - P0))]
                          + [wall_front(x, P[:, 0], P0 + f * (P1 - P0)) for f in WALL_LEVELS])
            if tp < T0 + 10e-6:
                continue
            xi = crossing(x, P[:, -1], 0.5 * (P0 + P1))
            T = triple_point(x, y, P, beta, xi, levels=(1.15, 1.25, 1.35))
            s, n = to_interface(np.array([p[0] for p in T]), np.array([p[1] for p in T]), beta)
            traj.append((tp, s[1], n[1], float(np.ptp(n))))
        fronts = np.array(fronts)
        m = np.isfinite(fronts[:, 1])
        r["front_speed_2mm"] = float(fit_line(fronts[m, 0], fronts[m, 1])[0])
        r["R_speed"] = U_S / np.sin(np.radians(beta))
        wall_speed = []
        for k in range(len(WALL_LEVELS)):
            xk = fronts[:, 2 + k]
            m = np.isfinite(xk) & (xk > 0.010) & (xk < 0.070)
            wall_speed.append(float(fit_line(fronts[m, 0], xk[m])[0]) if m.sum() > 2 else np.nan)
        r["wall_levels"] = {"f": list(WALL_LEVELS), "t": fronts[:, 0].tolist(),
                            "x": fronts[:, 2:].T.tolist(), "speed": wall_speed}
        r["S_t_speed_1pct"] = wall_speed[0]
        traj = np.array(traj)
        r["triple_point"] = {"t": traj[:, 0].tolist(), "s": traj[:, 1].tolist(), "n": traj[:, 2].tolist(),
                             "n_spread": traj[:, 3].tolist()}
        late = traj[:, 0] >= 0.5 * traj[-1, 0]
        b, e, _ = fit_line(traj[late, 1], traj[late, 2])
        r["chi_apex_last"] = float(np.degrees(np.arctan2(traj[-1, 2], traj[-1, 1])))
        r["chi_fit"], r["chi_fit_err"] = float(np.degrees(np.arctan(b))), float(np.degrees(e))
        r["n_T_last"] = float(traj[-1, 2])
        res["runs"][tag] = r

        line = (f"beta = {tag:>4}: U_s = {r['U_s']:.1f} m/s, transmitted front 2 mm deep {r['front_speed_2mm']:.0f} m/s "
                f"(R {r['R_speed']:.0f} m/s)")
        if "S_t_speed" in r:
            line += f", S_t speed = {r['S_t_speed']:.0f} +/- {r['S_t_err']:.0f} m/s"
        line += " (wall levels " + ", ".join(f"{100 * f:g}%: {v:.0f}" for f, v in zip(WALL_LEVELS, wall_speed)) + " m/s)"
        if "omega_r" in r:
            line += f", omega_r = {r['omega_r']:.2f} +/- {r['omega_r_err']:.2f} deg"
        if "phi_b" in r:
            line += f", phi_b = {r['phi_b']:.2f} +/- {r['phi_b_err']:.2f} deg"
        line += (f", T: n = {1e3 * r['n_T_last']:.2f} mm, chi (apex) = {r['chi_apex_last']:.2f} deg, "
                 f"chi (fit) = {r['chi_fit']:.2f} +/- {r['chi_fit_err']:.2f} deg")
        print(line)

    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print("wrote", out)
    return res


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def waves_figure(M, out):
    """Wave angles, precursor and triple point against polar theory and the experiment."""
    warnings.simplefilter("ignore", RuntimeWarning)
    col = {"15": "#4477AA", "20.5": "#228833", "30": "#E1701A", "50": "#6F4E9C"}
    gray = "0.55"
    b_rrr, b_det = np.degrees(polar.rrr_limit()), np.degrees(polar.detachment_angle())
    beta = np.linspace(2.0, b_det, 200)
    om = np.array([np.degrees(polar.regular_reflection(np.radians(b))[0]) for b in beta])
    beta_t = np.linspace(4.0, b_rrr - 1e-3, 120)
    phi = np.array([np.degrees(polar.rrr(np.radians(b))["phi_b"]) for b in beta_t])

    pp.use_style()
    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(10.0, 3.3), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.15, top=0.97, wspace=0.3)

    # (a) wave angles
    ax.plot(beta, om, color="0.35", lw=1.2, label=r"polar theory, $\omega_r$")
    ax.plot(beta_t, phi, color="0.35", lw=1.2, ls="--", label=r"polar theory, $\phi_b$")
    for b in (b_rrr, b_det):
        ax.axvline(b, color=gray, lw=0.8, zorder=0)
    kw = dict(ms=6, mew=1.1, lw=1.0, capsize=2.5, ls="none", zorder=5)
    for k, r in M.items():
        if "omega_r" in r:
            ax.errorbar(float(k), r["omega_r"], yerr=r["omega_r_err"], fmt="o", color=col[k], mfc="white", **kw)
        if "phi_b" in r:
            ax.errorbar(float(k), r["phi_b"], yerr=r["phi_b_err"], fmt="s", color=col[k], mfc="white", **kw)
    ax.errorbar([], [], fmt="o", color="0.2", mfc="white", ms=6, label=r"COMPAS, $\omega_r$ (reflected)")
    ax.errorbar([], [], fmt="s", color="0.2", mfc="white", ms=6, label=r"COMPAS, $\phi_b$ (transmitted)")
    ax.set_xlim(0, 55)
    ax.set_ylim(0, 100)
    for x, lab in ((0.5 * b_rrr, "RRR"), (0.5 * (b_rrr + b_det), "FPR"), (0.5 * (b_det + 55), "FMR")):
        ax.text(x, 96, lab, ha="center", va="top", fontsize=9, color=pp.INK2)
    ax.text(b_rrr + 0.6, 56, rf"${b_rrr:.2f}^\circ$", fontsize=8, color=pp.INK2, rotation=90, va="bottom")
    ax.text(b_det + 0.6, 3, rf"${b_det:.2f}^\circ$", fontsize=8, color=pp.INK2, rotation=90, va="bottom")
    ax.set_xlabel(r"$\beta$ (deg)")
    ax.set_ylabel("wave angle to the interface (deg)")
    ax.legend(loc="center right", fontsize=7.5, bbox_to_anchor=(1.0, 0.6))
    ax.text(0.97, 0.04, "(a)", transform=ax.transAxes, ha="right")

    # (b) leading edge of the precursor along the bottom wall
    bx.set_xlim(0, 60)
    bx.set_ylim(0, 85)
    bx.plot([0.0, 60.0], [0.0, 60e-3 * C_WATER], color=gray, lw=1.0, zorder=1)
    fig.canvas.draw()
    p0, p1 = bx.transData.transform([(0.0, 0.0), (60.0, 60e-3 * C_WATER)])
    bx.text(53.0, 1e-3 * C_WATER * 53.0 - 3.5, rf"$c_w (t - t_0)$, $c_w = {C_WATER:.0f}$ m/s", fontsize=7.5,
            color=pp.INK2, ha="right", va="top", rotation=np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0])),
            rotation_mode="anchor")
    sizes = {"15": 8.0, "20.5": 6.2, "30": 4.4, "50": 2.8}
    for k in M:
        W = M[k]["wall_levels"]
        t, x = 1e6 * (np.array(W["t"]) - T0), 1e3 * np.array(W["x"])
        m = np.isfinite(x[0]) & (x[0] < 79.0)
        bx.plot(t[m], x[0][m], "o", color=col[k], mfc="white" if k != "50" else col[k], ms=sizes[k], mew=1.0,
                zorder=5, label=rf"COMPAS, $\beta = {k}^\circ$")
    bx.set_xlabel(r"$t - t_0$ ($\mu$s)")
    bx.set_ylabel(r"precursor along the wall, $x$ (mm)")
    bx.legend(loc="upper left", fontsize=7.5)
    bx.text(0.97, 0.05, "(b)", transform=bx.transAxes, ha="right")

    # (c) triple point of the beta = 50 deg run
    T = M["50"]["triple_point"]
    s, n = 1e3 * np.array(T["s"]), 1e3 * np.array(T["n"])
    late = np.array(T["t"]) >= 0.5 * T["t"][-1]
    b, a = np.polyfit(s[late], n[late], 1)
    sf = np.array([s[late].min(), 55.0])              # the fitted stretch of the path
    n0 = a + b * sf[0]
    for c in (1.9, 2.1):
        cx.plot(sf, n0 + (sf - sf[0]) * np.tan(np.radians(c)), color="0.8", lw=0.8, zorder=0)
    h_exp, = cx.plot(sf, n0 + (sf - sf[0]) * np.tan(np.radians(2.0)), color=gray, lw=1.0, zorder=1,
                     label=r"$\chi = 2 \pm 0.1^\circ$, experiment")
    h_pts = cx.errorbar(s, n, yerr=0.5 * 1e3 * np.array(T["n_spread"]), fmt="o", color=col["50"], mfc="white",
                        ms=5, mew=1.0, lw=0.9, capsize=2, zorder=5, label=r"COMPAS, $\beta = 50^\circ$ (FMR)")
    h_fit, = cx.plot(sf, a + b * sf, color=col["50"], lw=0.9,
                     label=rf"COMPAS fit, $\chi = {np.degrees(np.arctan(b)):.1f}^\circ$")
    cx.set_xlim(0, 55)
    cx.set_ylim(0, 0.5 * np.ceil(2.5 * n.max()))
    cx.set_xlabel(r"$s_T$ along the interface (mm)")
    cx.set_ylabel(r"$n_T$ above the interface (mm)")
    cx.legend(handles=[h_pts, h_fit, h_exp], loc="upper left", fontsize=7.5)
    cx.text(0.97, 0.05, "(c)", transform=cx.transAxes, ha="right")
    pp.save(fig, out, tight=False)


def snapshots_figure(root, out):
    """The four refraction patterns at the last plotfile of each run."""
    g_rho, g_p, darkest = (60.0, 6000.0), (1.0e7, 4.0e8), 0.45    # schlieren limits, kg/m^4 and Pa/m
    pp.use_style()
    FW = 10.0
    fig, axes = pp.panel_grid(2, 2, 0.4, FW, left=0.62, right=0.95, top=0.1, bottom=0.52, wgap=0.14, hgap=0.14)
    label = dict(transform=None, va="top", fontsize=10, color=pp.INK, zorder=10,
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    hs = []                                             # legend handles of the AMR levels
    for k, (beta, regime) in enumerate(RUNS):
        ax = axes.flat[k]
        fr = pp.load(pp.pick(run_dir(root, beta))[1], ("a1", "rho", "P"))
        A = np.clip(fr.q["a1"], 0, 1)
        air = minimum_filter(A, size=5) > 0.99        # pure air, also in the cells next to it
        s = np.where(air, pp.schlieren(fr.q["rho"], fr.dx, lims=g_rho, darkest=darkest),
                     pp.schlieren(fr.q["P"], fr.dx, lims=g_p, darkest=darkest))
        img = pp.TEAL_SAND_LIGHT(A.T)[..., :3] * (1.0 - s.T[..., None])
        ext = 1e3 * np.array(fr.extent)
        ax.imshow(img, origin="lower", extent=ext, interpolation="antialiased")
        ax.contour(1e3 * fr.x, 1e3 * fr.y, A.T, levels=[0.5], colors=pp.INK, linewidths=0.35)
        if beta == 50.0:
            hs = pp.draw_boxes(ax, fr, 1e3, lw=0.5, alpha=0.9, zorder=4)
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_aspect("equal")
        label["transform"] = ax.transAxes
        ax.text(0.012, 0.965, rf"$\beta = {beta:g}^\circ$, {regime}", ha="left", **label)
        ax.text(0.988, 0.965, rf"$t - t_0 = {1e6 * (fr.t - T0):.1f}\ \mu$s", ha="right", **label)
        ax.set_xticks([0, 20, 40, 60])
        ax.set_yticks([0, 20, 40] if k < 2 else [0, 20])
        if k >= 2:
            ax.set_xlabel(r"$x$ (mm)")
        else:
            ax.tick_params(labelbottom=False)
        if k % 2 == 0:
            ax.set_ylabel(r"$y$ (mm)")
        else:
            ax.tick_params(labelleft=False)
    t_, b_ = axes[0, 1].get_position().y1, axes[0, 1].get_position().y0     # color bar beside the upper row
    cax = fig.add_axes([axes[0, 1].get_position().x1 + 0.15 / FW, b_, 0.12 / FW, t_ - b_])
    cb = fig.colorbar(plt.cm.ScalarMappable(Normalize(0, 1), pp.TEAL_SAND_LIGHT), cax=cax, ticks=[0, 0.5, 1])
    cb.set_ticklabels(["0 water", "0.5", "1 air"])
    cb.set_label(r"$\alpha_1$")
    cb.ax.minorticks_off()
    cb.ax.axhline(0.5, color=pp.INK, lw=0.35)
    if hs:                                              # AMR levels, beside the beta = 50 deg panel
        fig.legend(hs, [f"level {L + 1}" for L in range(len(hs))], loc="upper left", handlelength=1.4,
                   bbox_to_anchor=(axes[1, 1].get_position().x1 + 0.08 / FW, axes[1, 1].get_position().y1),
                   borderaxespad=0, labelspacing=0.35)
    pp.save(fig, out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["waves", "snapshots"])
    ap.add_argument("root", help="directory with one run per angle, b<beta>/plot")
    ap.add_argument("--json", help="measurements file written by waves (default <root>/measurements.json)")
    ap.add_argument("--out", help="figure file")
    args = ap.parse_args()
    if args.command == "waves":
        js = Path(args.json) if args.json else Path(args.root) / "measurements.json"
        res = measure(args.root, [b for b, _ in RUNS], js)
        waves_figure(res["runs"], args.out or "shockrefraction_waves.png")
    else:
        snapshots_figure(args.root, args.out or "shockrefraction_snapshots.png")
