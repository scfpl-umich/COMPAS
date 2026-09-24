# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Shock polars for a stiffened gas, and the regular refraction of a shock at an air-water
interface, after Appendix A and section 3 of Anbu Serene Raj et al. (J. Fluid Mech. 998, A49, 2024).

In the frame of the refraction point R, which runs up the interface at U_s / sin(beta), the air and
the water ahead of the waves both flow along the interface, the air at M0 = Ms / sin(beta) and the
water at Mb = M0 a0 / ab. The incident shock makes the angle beta with this flow. An oblique shock
of angle phi in a stiffened gas deflects the flow by the same theta(M, phi) as in an ideal gas
(eq. A17), and raises the pressure by

    p2/p1 = 1 + 2 gamma/(gamma + 1) (Mn^2 - 1) (1 + pinf/p1),   Mn = M sin(phi)      (eq. A15)

where the Mach numbers use the stiffened sound speed sqrt(gamma (p + pinf)/rho). In regular
refraction with a transmitted shock (RRR), the incident and reflected shocks turn the air by
theta1 - theta2, the transmitted shock turns the water by the same angle, and the pressures behind
the reflected and the transmitted shocks are equal. RRR ends when the reflected polar is tangent to
the transmitted polar at its normal-shock point (eqs. 3.5 and 3.6), just before Mb = 1, where the
trace speed U_s / sin(beta) equals the sound speed of the water. Beyond it the transmitted wave is
a free precursor, and the reflection in the air is regular (FPR) until the reflected shock can no
longer turn the flow back (detachment), and a Mach reflection (FMR) after that.

    python3 analysis/polar.py           prints the transition angles and the wave angles
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import postprocess                                      # noqa: E402,F401  (sets up the Python environment)

import numpy as np                                      # noqa: E402
from scipy.optimize import brentq                       # noqa: E402

# Air (ideal gas) and water (stiffened gas) at rest, as in the case inputs (SI units)
AIR   = dict(gamma=1.4, pinf=0.0,   rho=1.17665)
WATER = dict(gamma=2.8, pinf=8.5e8, rho=1053.018)
P0    = 101325.0
MS    = 1.46


def sound_speed(gas, p=P0):
    return np.sqrt(gas["gamma"] * (p + gas["pinf"]) / gas["rho"])


def oblique(M, phi, gas, p1):
    """Deflection theta (rad), pressure ratio p2/p1 and downstream Mach number M2 of an oblique
    shock at angle phi (rad) to a flow of Mach number M at pressure p1 (eqs. A10, A15, A17)."""
    g = gas["gamma"]
    Mn2 = (M * np.sin(phi))**2
    theta = np.arctan(2.0 / np.tan(phi) * (Mn2 - 1.0) / (M**2 * (g + np.cos(2.0 * phi)) + 2.0))
    pr = 1.0 + 2.0 * g / (g + 1.0) * (Mn2 - 1.0) * (1.0 + gas["pinf"] / p1)
    M2n = np.sqrt((1.0 + 0.5 * (g - 1.0) * Mn2) / (g * Mn2 - 0.5 * (g - 1.0)))
    return theta, pr, M2n / np.sin(phi - theta)


def max_deflection(M, gas, p1):
    """Shock angle and deflection at the detachment point of the polar of Mach number M."""
    phi = np.linspace(np.arcsin(1.0 / M) + 1e-9, 0.5 * np.pi - 1e-9, 20001)
    theta = oblique(M, phi, gas, p1)[0]
    k = np.argmax(theta)
    return phi[k], theta[k]


def weak_angle(M, theta, gas, p1):
    """Shock angle of the weak solution that deflects a flow of Mach number M by theta >= 0,
    or nan beyond detachment."""
    phi_d, theta_d = max_deflection(M, gas, p1)
    if theta > theta_d:
        return np.nan
    mu = np.arcsin(1.0 / M)
    if theta <= 0.0:
        return mu
    return brentq(lambda p: oblique(M, p, gas, p1)[0] - theta, mu + 1e-12, phi_d)


def incident(beta, Ms=MS, air=AIR, water=WATER):
    """Mach numbers M0, Mb in the frame of R, and the deflection theta1, pressure p1 and Mach
    number M1 behind the incident shock, for the interface angle beta (rad)."""
    M0 = Ms / np.sin(beta)
    Mb = M0 * sound_speed(air) / sound_speed(water)
    theta1, pr1, M1 = oblique(M0, beta, air, P0)
    return M0, Mb, theta1, pr1 * P0, M1


def regular_reflection(beta, Ms=MS, air=AIR, water=WATER):
    """Reflection in the air as on a solid wall (theta2 = theta1): the reflected shock angle
    omega_r = phi2 - theta2 to the interface (rad), and the pressure p2 behind it.
    nan beyond detachment."""
    M0, Mb, theta1, p1, M1 = incident(beta, Ms, air, water)
    phi2 = weak_angle(M1, theta1, air, p1)
    if np.isnan(phi2):
        return np.nan, np.nan
    return phi2 - theta1, p1 * oblique(M1, phi2, air, p1)[1]


def rrr(beta, Ms=MS, air=AIR, water=WATER):
    """Regular refraction with a transmitted shock at the interface angle beta (rad): the
    transmitted shock angle phi_b, the reflected shock angle omega_r = phi2 - theta2 (both to the
    interface, rad), the interface deflection delta (rad) and the pressure behind the waves.
    Returns None when there is no RRR solution."""
    M0, Mb, theta1, p1, M1 = incident(beta, Ms, air, water)
    if Mb <= 1.0:
        return None
    mu_b = np.arcsin(1.0 / Mb)

    def mismatch(phi_b):
        # reflected-shock pressure minus transmitted-shock pressure at the same interface deflection
        delta, prb, _ = oblique(Mb, phi_b, water, P0)
        phi2 = weak_angle(M1, theta1 - delta, air, p1)
        return p1 * oblique(M1, phi2, air, p1)[1] - prb * P0

    # The transmitted pressure rises with phi_b along the whole t-polar, up to the normal shock;
    # close to the limit the solution lies past the detachment point of the t-polar
    if mismatch(0.5 * np.pi) > 0.0:
        return None
    phi_b = brentq(mismatch, mu_b + 1e-12, 0.5 * np.pi, xtol=1e-14)
    delta, prb, _ = oblique(Mb, phi_b, water, P0)
    phi2 = weak_angle(M1, theta1 - delta, air, p1)
    return dict(phi_b=phi_b, omega_r=phi2 - (theta1 - delta), delta=delta, p=prb * P0)


def rrr_limit(Ms=MS, air=AIR, water=WATER):
    """Interface angle (rad) of the RRR to BPR transition: the reflected polar at zero net
    deflection reaches the normal-shock pressure of the transmitted polar (eqs. 3.5, 3.6)."""
    def f(beta):
        M0, Mb, theta1, p1, M1 = incident(beta, Ms, air, water)
        p2 = regular_reflection(beta, Ms, air, water)[1]
        p_ns = P0 * oblique(Mb, 0.5 * np.pi, water, P0)[1]
        return p_ns - p2
    b_sonic = sonic_angle(Ms, air, water)
    return brentq(f, np.radians(5.0), b_sonic - 1e-12, xtol=1e-15)


def sonic_angle(Ms=MS, air=AIR, water=WATER):
    """Interface angle (rad) at which Mb = 1, i.e. U_s / sin(beta) = ab."""
    return np.arcsin(Ms * sound_speed(air) / sound_speed(water))


def detachment_angle(Ms=MS, air=AIR, water=WATER):
    """Interface angle (rad) beyond which the reflected shock cannot turn the air back parallel to
    the interface, the detachment criterion for regular reflection (FPR to FMR)."""
    def f(beta):
        M0, Mb, theta1, p1, M1 = incident(beta, Ms, air, water)
        return max_deflection(M1, air, p1)[1] - theta1
    return brentq(f, np.radians(25.0), np.radians(55.0), xtol=1e-12)


if __name__ == "__main__":
    deg = np.degrees
    print(f"air sound speed {sound_speed(AIR):.2f} m/s, water sound speed {sound_speed(WATER):.2f} m/s, "
          f"U_s = {MS * sound_speed(AIR):.2f} m/s")
    print(f"RRR -> BPR (polar tangency): beta = {deg(rrr_limit()):.4f} deg; Mb = 1 at beta = {deg(sonic_angle()):.4f} deg")
    print(f"FPR -> FMR (detachment in the air): beta = {deg(detachment_angle()):.3f} deg")

    # The same with the paper's Table 1 air density and the Mach number of its Table 2 state
    # (p1 = 235173.8 Pa over 101325 Pa, Ms = 1.46023), for comparison with its 19.7148 deg
    air_t1 = dict(AIR, rho=1.176)
    Ms_t2 = np.sqrt(1.0 + (235173.823 / P0 - 1.0) * 2.4 / 2.8)
    print(f"with rho_air = 1.176 and Ms = {Ms_t2:.5f}: RRR -> BPR at {deg(rrr_limit(Ms_t2, air_t1)):.4f} deg, "
          f"Mb = 1 at {deg(sonic_angle(Ms_t2, air_t1)):.4f} deg")

    for b in (15.0, 20.5, 30.0, 50.0):
        M0, Mb, theta1, p1, M1 = incident(np.radians(b))
        s = rrr(np.radians(b))
        om, p2 = regular_reflection(np.radians(b))
        line = f"beta = {b:5.1f}: M0 = {M0:.3f}, M1 = {M1:.3f}, Mb = {Mb:.3f}"
        if s is not None:
            line += (f"  RRR: phi_b = {deg(s['phi_b']):.3f} deg, omega_r = {deg(s['omega_r']):.3f} deg, "
                     f"delta = {deg(s['delta']):.4f} deg, p/p0 = {s['p'] / P0:.3f}")
        elif not np.isnan(om):
            line += f"  FPR: omega_r = {deg(om):.3f} deg, p2/p0 = {p2 / P0:.3f}"
        else:
            line += "  FMR (no regular reflection)"
        print(line)
