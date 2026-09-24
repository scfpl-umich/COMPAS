# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""Contour plots of the Guderley implosion.

The simulated quarter plane (x, y >= 0) is mirrored about both axes. The upper half of each
panel shows the density, and the lower half shows the pressure with the interface
(alpha_1 = 0.5) drawn in white. The color ranges stop at the 99.5th percentile, so that the
peaks at the focusing point do not wash out the rest of the panel.

    python analysis/plot_contours.py plot/Guderley --times 0.1 0.3 0.44 0.52 0.7 0.9
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compas_io import find_plotfiles, plotfile_time, load_uniform

parser = argparse.ArgumentParser()
parser.add_argument("plot_dir", help="directory with the plotfiles, e.g. plot/Guderley")
parser.add_argument("--times", type=float, nargs="+", default=[0.1, 0.3, 0.44, 0.52, 0.7, 0.9])
parser.add_argument("--extent", type=float, default=0.62, help="half-width of the plotted region")
parser.add_argument("--out", default="contours.png")
args = parser.parse_args()

# Pick the plotfile closest to each requested time
plotfiles = find_plotfiles(args.plot_dir)
times = np.array([plotfile_time(p) for p in plotfiles])
chosen = [plotfiles[int(np.argmin(np.abs(times - t)))] for t in args.times]

ncol = 3
nrow = int(np.ceil(len(chosen) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.7 * ncol, 3.1 * nrow), constrained_layout=True)

for ax, plotfile in zip(np.ravel(axes), chosen):
    t, x, y, d = load_uniform(plotfile, fields=("a1", "rho", "P"))
    L = x[-1] + 0.5 * (x[1] - x[0])

    # Upper half: density, mirrored across the y axis (arrays are indexed [i, j], i along x)
    rho = np.concatenate([d["rho"][::-1, :], d["rho"]], axis=0)
    im_rho = ax.imshow(rho.T, origin="lower", extent=[-L, L, 0, L], cmap="viridis",
                       vmin=rho.min(), vmax=np.percentile(rho, 99.5))

    # Lower half: pressure, mirrored across both axes, with the interface
    P = np.concatenate([d["P"][::-1, ::-1], d["P"][:, ::-1]], axis=0)
    im_P = ax.imshow(P.T, origin="lower", extent=[-L, L, -L, 0], cmap="magma",
                     vmin=P.min(), vmax=np.percentile(P, 99.5))
    a1 = np.concatenate([d["a1"][::-1, ::-1], d["a1"][:, ::-1]], axis=0)
    xf = np.concatenate([-x[::-1], x])
    ax.contour(xf, -y[::-1], a1.T, levels=[0.5], colors="w", linewidths=0.6)

    # One colorbar beside each half
    for im, box, label in [(im_rho, [1.03, 0.52, 0.05, 0.46], r"$\rho$"),
                           (im_P,   [1.03, 0.02, 0.05, 0.46], r"$P$")]:
        cb = fig.colorbar(im, cax=ax.inset_axes(box), extend="max")
        cb.set_label(label, rotation=0, labelpad=8)
        cb.ax.tick_params(labelsize=7)

    ax.set_title(f"t = {t:.2f}")
    ax.set_xlim(-args.extent, args.extent)
    ax.set_ylim(-args.extent, args.extent)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

for ax in np.ravel(axes)[len(chosen):]:
    ax.axis("off")

fig.savefig(args.out, dpi=200)
print(f"Wrote {args.out}")
