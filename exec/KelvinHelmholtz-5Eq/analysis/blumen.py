# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Linear theory of the compressible tanh shear layer (Blumen 1970).

Source: W. Blumen, Shear layer instability of an inviscid compressible fluid, J. Fluid Mech. 40
(1970) 769-781. The basic flow is u = U tanh(y/L) at uniform pressure and density, lengths are in
units of L, velocities in units of U and times in units of L/U, and the Mach number is M = U/a,
with a the (uniform) sound speed. Disturbances go as exp[i alpha (x - c t)], so the growth rate is
alpha c_i. TABLE1 is Blumen's Table 1, the fastest-growing mode (c_r = 0) at each M: its
wavenumber alpha and its growth rate alpha c_i.

The solver integrates Blumen's pressure equation (7),
    (u - c) p'' - 2 u' p' - alpha^2 (u - c) [1 - M^2 (u - c)^2] p = 0,
for c = i c_i. With c_r = 0 the real part of p is even in y and the imaginary part odd
(Blumen, section 5), so only y >= 0 is needed: the Riccati variable G = p'/p is integrated from
the decaying solution at large y (or p' = 0 at a slip wall, Blumen's eq. 12) down to y = 0,
where symmetry requires Re G(0) = 0. This checks Table 1 and gives the rate for a finite
channel of half-height Y.

    python3 analysis/blumen.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess                                      # noqa: E402,F401  (sets up the Python environment)

import numpy as np                                      # noqa: E402
from scipy.integrate import solve_ivp                   # noqa: E402
from scipy.optimize import brentq, minimize_scalar      # noqa: E402

# Blumen (1970), Table 1: M, alpha of the fastest-growing mode, and its growth rate alpha c_i
TABLE1 = np.array([
    # M     alpha   alpha c_i
    [0.0,  0.445,  0.190],
    [0.1,  0.433,  0.187],
    [0.2,  0.426,  0.181],
    [0.3,  0.417,  0.171],
    [0.4,  0.409,  0.158],
    [0.5,  0.397,  0.141],
    [0.6,  0.370,  0.122],
    [0.7,  0.326,  0.101],
    [0.8,  0.279,  0.078],
    [0.9,  0.208,  0.055],
    [1.0,  0.000,  0.000],
])


def mismatch(ci, alpha, M, Y=None, y_far=30.0):
    """Re G(0) for c = i ci, integrating G = p'/p from y_far (decaying solution) or from a
    slip wall at y = Y (p' = 0) down to y = 0."""
    c = 1j * ci
    if Y is None:
        y0 = y_far
        G0 = -alpha * np.sqrt(1.0 - M**2 * (1.0 - c)**2 + 0j)   # decaying branch
        if G0.real > 0.0:
            G0 = -G0
    else:
        y0, G0 = Y, 0.0 + 0.0j

    def rhs(y, G):
        W = np.tanh(y) - c
        du = 1.0 / np.cosh(y)**2
        return 2.0 * du / W * G + alpha**2 * (1.0 - M**2 * W**2) - G**2

    sol = solve_ivp(rhs, (y0, 0.0), [complex(G0)], method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.y[0, -1].real


def growth_rate(alpha, M, Y=None):
    """Growth rate alpha c_i of the c_r = 0 mode at wavenumber alpha (0 if it is stable)."""
    if alpha**2 + M**2 >= 1.0:
        return 0.0
    grid = np.linspace(1e-3, 1.0, 60)
    f = [mismatch(ci, alpha, M, Y) for ci in grid]
    for a, b, fa, fb in zip(grid[::-1][1:], grid[::-1][:-1], f[::-1][1:], f[::-1][:-1]):
        if fa * fb < 0.0:                      # the largest root is the unstable mode
            return alpha * brentq(mismatch, a, b, args=(alpha, M, Y), xtol=1e-12)
    return 0.0


def most_unstable(M, Y=None):
    """(alpha_max, growth rate) of the fastest-growing c_r = 0 mode at Mach number M."""
    hi = np.sqrt(max(1.0 - M**2, 1e-6))
    r = minimize_scalar(lambda a: -growth_rate(a, M, Y), bounds=(0.05 * hi, 0.95 * hi),
                        method="bounded", options={"xatol": 1e-5})
    return r.x, -r.fun


if __name__ == "__main__":
    print(" M    alpha_max (Blumen, here)   alpha c_i (Blumen, here)   rate at Blumen's alpha, walls at |y| = 20, 40")
    for M, a_b, s_b in TABLE1[:-1]:
        a, s = most_unstable(M)
        s20 = growth_rate(a_b, M, Y=20.0)
        s40 = growth_rate(a_b, M, Y=40.0)
        print(f"{M:.1f}   {a_b:.3f}  {a:.4f}              {s_b:.3f}  {s:.4f}"
              f"              {growth_rate(a_b, M):.4f}  {s20:.4f}  {s40:.4f}")
