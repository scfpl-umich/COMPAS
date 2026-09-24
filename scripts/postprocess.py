#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
# Copyright (C) 2026 The Regents of the University of Michigan
# This file is part of COMPAS, free software under the GNU General Public License,
# version 3 or later, with an additional permission under section 7. It comes with
# ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

"""
===============================================================================
COMPAS post-processing: plotfiles, snapshot figures and movies
===============================================================================

Reads the AMReX plotfiles of 2D COMPAS runs with yt and draws them in the style of the COMPAS
documentation. It is a command-line tool and a module that the analysis scripts of the cases
import.

A run is a plot directory, the directory with the plotfiles (plt*), for example
plot/KH-M0.4. The copies AMReX leaves when it overwrites a plotfile (*.old.*) are skipped.

-------------------------------------------------------------------------------
COMMANDS
-------------------------------------------------------------------------------

list RUN [RUN ...]
    Print the plotfiles of each run and their times.

snapshots RUN [RUN ...] [--times T ...]
    One row of panels per run and one column per time (the plotfile nearest each time, or the
    last plotfile without --times). Writes --out (about 4800 px wide) and a 400 dpi copy with
    the suffix _400dpi.

movie RUN [RUN ...] [--tmin T] [--tmax T] [--every N] [--sync time|span] [--fps N]
    The same panels over every plotfile of the first run between tmin and tmax, one row per run,
    encoded by ffmpeg as H.264 (yuv420p). Each run is at its plotfile nearest the frame time, or
    with --sync span at the same fraction of the time span of its own plotfiles in [tmin, tmax],
    for runs that end at different times. Writes --out (default movie.mp4, 1600 px wide) and the
    last frame as <out>_poster.png.

interface RUN [RUN ...] [--times T ...] [--field a1] [--xlim X0 X1] [--ylim Y0 Y1]
    Measures of a volume-fraction interface in each plotfile: the extremes of the field, the
    length of its 0.5 contour, and the 5 to 95 percent thickness, the area where
    0.05 < field < 0.95 over that length, in cells of the finest level.

Drawing options of snapshots and movie:

    --field NAME            field drawn in color (default a1; none for a white background)
    --cmap NAME             teal-sand (default), teal-sand-light, or a matplotlib colormap
    --range LO HI           color range (default 0 1 for a1, else the range over the panels)
    --schlieren FIELD       darken the colors with a light schlieren of |grad FIELD|
    --schlieren-range G0 G1 |grad| where the shading starts and where it is darkest, on a log
                            scale in between (default from the panels, printed)
    --contour FIELD LEVELS  contour lines of FIELD at the given values; START:STOP:N gives N
                            equally spaced values
    --amr [RUN ...]         outline the AMR levels, in every panel or in the panels of the
                            given runs (counted from 1), in the color of each level
    --amr-style regions|boxes
                            with --amr, the outline of the region each level covers (default),
                            or every patch box of each level
    --amr-lw W, --amr-alpha A
                            with --amr, the width in points and the opacity of the AMR lines
                            (default 0.4 and 0.5)
    --xlim X0 X1, --ylim Y0 Y1   the window (default the domain)
    --mirror x|y            reflect each panel about the lower x or y edge of the window, a
                            symmetry plane of the domain
    --split FIELD           with --mirror, draw FIELD in the reflected half (--split-cmap,
                            --split-range), with its own color bar
    --scale S               multiply lengths by S on the axes (for example 1e3 for mm)
    --t0 T, --tscale S      show times as (t - T) S
    --title FMT             title of each column, e.g. '$t = {t:.0f}$'
    --label FMT [FMT ...]   label in each panel, one for all runs or one per run; the default,
                            when no title is given, is '$t = {t:.3g}$'
    --row-label TEXT ...    label to the left of each row
    --xlabel, --ylabel      axis labels (default $x$, $y$)
    --colorbar TEXT         label of the color bar (default the field name)
    --caption LINE ...      caption above the panels, one argument per line
    --wrap N                with one time, the runs in rows of N panels instead of one per row
    --interpolation NAME    matplotlib interpolation of the images (default antialiased, which
                            shows each cell as a block once it spans several pixels; bilinear
                            blends them)
    --frame light|box       light gray frame with outward ticks (default), or the boxed axes
                            of the line plots
    --width W               figure width in inches (default 10)
    --width-px N            width of the web figure or movie in pixels

Examples, from a case directory:

    python3 ../../scripts/postprocess.py list plot/KH-M0.4
    python3 ../../scripts/postprocess.py snapshots plot/KH-M0.4 --times 65 75 85 95 --ylim -8.75 8.75
    python3 ../../scripts/postprocess.py movie plot/KH-M0.4 plot/KH-M0.4-THINC --tmin 40 --fps 10
    python3 ../../scripts/postprocess.py interface plot/KH-M0.4 --times 65 95

-------------------------------------------------------------------------------
MODULE
-------------------------------------------------------------------------------

The analysis scripts of the cases put scripts/ on sys.path and import this file:

    plotfiles(run)            [(time, path)], sorted by time
    pick(run, t)              (time, path) of the plotfile nearest t
    load(plotfile, fields, xlim, ylim)
                              Frame: the fields composited on the grid of the finest level,
                              each level upsampled linearly and used where no finer level
                              covers it, with the finest level present in each cell and the
                              patch boxes
    teal_sand(start)          the volume-fraction colormap, from deep teal (0) to light sand (1);
                              TEAL_SAND and its lighter part TEAL_SAND_LIGHT
    schlieren(q, dx, dy, lims, darkest)
                              shading in [0, darkest] from |grad q|, logarithmic between lims
    level_outlines(frame, level), draw_levels(ax, frame)
                              the outline of the region each AMR level covers
    box_outlines(frame, level), draw_boxes(ax, frame)
                              the outline of every patch box of each AMR level
    interface_measures(frame, field)
                              extremes, 0.5-contour length and 5 to 95 percent thickness
    panel_grid(...), light_frame(ax), save(fig, out), Movie(out)
                              layout, frame style, figure output and movie output

-------------------------------------------------------------------------------
PYTHON ENVIRONMENT
-------------------------------------------------------------------------------

Run outside a virtual environment, this script, and any script that imports it, creates or
reuses scripts/compas_python_env, installs the packages of scripts/requirements.txt that are
missing, and runs again inside it. ffmpeg must be on the PATH for movies.
===============================================================================
"""

import argparse
import os
import subprocess
import sys
import venv
from dataclasses import dataclass, field as dc_field
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ENV_DIR = SCRIPTS / "compas_python_env"
REQUIREMENTS = SCRIPTS / "requirements.txt"
STYLE = SCRIPTS / "compas.mplstyle"


# ---------------------------------------------------------------------------
# Python environment
# ---------------------------------------------------------------------------

def running_inside_venv():
    """True if already running inside a virtual environment."""
    return sys.prefix != sys.base_prefix or bool(os.environ.get("VIRTUAL_ENV"))


def ensure_env():
    """Create or reuse scripts/compas_python_env with the packages of requirements.txt, and
    return its Python executable."""
    if not ENV_DIR.exists():
        print(f"[COMPAS] Creating the Python environment at '{ENV_DIR}'")
        venv.EnvBuilder(with_pip=True).create(ENV_DIR)
    python = ENV_DIR / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    packages = [ln.split("#")[0].strip() for ln in REQUIREMENTS.read_text().splitlines()]
    missing = [p for p in packages if p and subprocess.run(
        [str(python), "-c", f"import {p}"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode]
    if missing:
        print(f"[COMPAS] Installing {' '.join(missing)}")
        subprocess.check_call([str(python), "-m", "pip", "install", *missing])
    return python


def bootstrap():
    """Run the calling script again inside scripts/compas_python_env when it was started
    outside a virtual environment."""
    script = sys.argv[0] if sys.argv else ""
    if running_inside_venv() or not script or not Path(script).is_file():
        return
    python = ensure_env()
    env = dict(os.environ, VIRTUAL_ENV=str(ENV_DIR), PATH=str(python.parent) + os.pathsep + os.environ.get("PATH", ""))
    sys.exit(subprocess.call([str(python), script, *sys.argv[1:]], env=env))


bootstrap()

import numpy as np                                            # noqa: E402
import matplotlib                                             # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                               # noqa: E402
from matplotlib.collections import LineCollection             # noqa: E402
from matplotlib.colors import ListedColormap, Normalize       # noqa: E402
from matplotlib.lines import Line2D                           # noqa: E402
from matplotlib.ticker import MaxNLocator                     # noqa: E402
from matplotlib.transforms import ScaledTranslation           # noqa: E402

# Colors of the figures: text, secondary text and the light panel frame; AMR levels (Okabe-Ito)
INK, INK2, FRAME = "#1f1e1c", "#6b6964", "#9a9892"
LEVEL_COLORS = {1: "#56B4E9", 2: "#E69F00", 3: "#D55E00", 4: "#CC79A7", 5: "#009E73"}


def use_style():
    """The matplotlib style of the COMPAS figures: serif Computer Modern, boxed axes."""
    plt.style.use(STYLE)


# ---------------------------------------------------------------------------
# Plotfiles
# ---------------------------------------------------------------------------

def plotfile_time(plotfile):
    """Simulation time of a plotfile, from its Header."""
    lines = (Path(plotfile) / "Header").read_text().splitlines()
    return float(lines[int(lines[1]) + 3])


def plotfiles(run, tmin=-np.inf, tmax=np.inf):
    """(time, path) of the plotfiles of a run (a plot directory, or one plotfile), sorted by
    time, between tmin and tmax."""
    run = Path(run)
    if (run / "Header").is_file():
        paths = [run]
    else:
        paths = [p for p in run.glob("plt*") if ".old" not in p.name and (p / "Header").is_file()]
    if not paths:
        raise SystemExit(f"no plotfiles in {run}")
    return [(t, p) for t, p in sorted((plotfile_time(p), p) for p in paths) if tmin <= t <= tmax]


def pick(run, t=None):
    """(time, path) of the plotfile of run nearest t, or the last one if t is None."""
    files = plotfiles(run)
    return files[-1] if t is None else min(files, key=lambda f: abs(f[0] - t))


@dataclass
class Frame:
    """Fields of one plotfile composited on the grid of its finest level. q[name] and level are
    indexed [i, j]; x and y are the cell centers, extent the edges [x0, x1, y0, y1]; boxes are
    the patches above level 0, as (level, x0, y0, x1, y1)."""
    t: float
    x: np.ndarray
    y: np.ndarray
    dx: float
    dy: float
    q: dict
    level: np.ndarray
    boxes: list = dc_field(default_factory=list)

    @property
    def extent(self):
        return [self.x[0] - self.dx / 2, self.x[-1] + self.dx / 2, self.y[0] - self.dy / 2, self.y[-1] + self.dy / 2]


def load(plotfile, fields=("a1",), xlim=None, ylim=None):
    """Load fields of a 2D plotfile as a Frame on the grid of its finest level, over the whole
    domain or over the window xlim, ylim widened to whole cells of level 0. Each level is
    upsampled linearly and used where no finer level covers it."""
    import yt
    from scipy.ndimage import zoom
    yt.set_log_level(40)
    ds = yt.load(str(plotfile))
    top, r = ds.max_level, ds.refine_by
    lo, hi = ds.domain_left_edge.d, ds.domain_right_edge.d
    n0 = ds.domain_dimensions[:2]
    d0 = (hi[:2] - lo[:2]) / n0

    def cells(lim, k):
        if lim is None:
            return 0, int(n0[k])
        a = int(np.floor((lim[0] - lo[k]) / d0[k] + 1e-9))
        b = int(np.ceil((lim[1] - lo[k]) / d0[k] - 1e-9))
        return max(a, 0), min(b, int(n0[k]))

    (i0, i1), (j0, j1) = cells(xlim, 0), cells(ylim, 1)
    f = r**top
    nx, ny = (i1 - i0) * f, (j1 - j0) * f
    dx, dy = d0 / f
    x0, y0 = lo[0] + i0 * d0[0], lo[1] + j0 * d0[1]

    level = np.zeros((nx, ny), np.int8)
    boxes = []
    for g in ds.index.grids:
        a = np.round((g.LeftEdge.d[:2] - (x0, y0)) / (dx, dy)).astype(int)
        b = np.round((g.RightEdge.d[:2] - (x0, y0)) / (dx, dy)).astype(int)
        if g.Level > 0:
            boxes.append((g.Level, *g.LeftEdge.d[:2], *g.RightEdge.d[:2]))
        s = level[max(a[0], 0):max(b[0], 0), max(a[1], 0):max(b[1], 0)]
        np.maximum(s, g.Level, out=s)

    q = {k: np.zeros((nx, ny)) for k in fields}
    for L in range(top + 1):
        m = level == L
        if not m.any():
            continue
        c = r**(top - L)
        cg = ds.covering_grid(level=L, left_edge=[x0, y0, lo[2]], dims=[nx // c, ny // c, 1])
        for k in fields:
            v = cg["boxlib", k][:, :, 0].d
            if c > 1:
                v = zoom(v, c, order=1, mode="nearest", grid_mode=True)
            q[k][m] = v[m]
    x = x0 + (np.arange(nx) + 0.5) * dx
    y = y0 + (np.arange(ny) + 0.5) * dy
    return Frame(float(ds.current_time), x, y, dx, dy, q, level, boxes)


# ---------------------------------------------------------------------------
# Colors, schlieren and AMR levels
# ---------------------------------------------------------------------------

_M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929], [0.2119034982, 0.6806995451, 0.1073969566],
                [0.0883024619, 0.2817188376, 0.6299787005]])
_M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468], [1.9779984951, -2.4285922050, 0.4505937099],
                [0.0259040371, 0.7827717662, -0.8086757660]])
# (L, C, h) anchors in OKLCH from deep teal (0) to light sand (1)
TEAL_SAND_ANCHORS = np.array([(0.30, 0.085, 252), (0.50, 0.080, 215), (0.72, 0.060, 150), (0.86, 0.070, 95),
                              (0.95, 0.045, 85)])


def oklch_to_rgb(L, C, h):
    """sRGB color of an OKLCH color."""
    lab = np.array([L, C * np.cos(np.radians(h)), C * np.sin(np.radians(h))])
    lin = np.linalg.inv(_M1) @ (np.linalg.inv(_M2) @ lab)**3
    return np.clip(np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * np.abs(lin)**(1 / 2.4) - 0.055), 0, 1)


def teal_sand(start=0.0, n=256):
    """Sequential teal-to-sand colormap whose lightness rises linearly, from the fraction start
    of the full map to its light end."""
    A = TEAL_SAND_ANCHORS
    Ls = np.linspace(A[0, 0] + start * (A[-1, 0] - A[0, 0]), A[-1, 0], n)
    return ListedColormap([oklch_to_rgb(l, np.interp(l, A[:, 0], A[:, 1]), np.interp(l, A[:, 0], A[:, 2]))
                           for l in Ls], name="teal_sand" if start == 0 else "teal_sand_light")


TEAL_SAND = teal_sand()
# Without the darkest 0.4 of the map, so that a schlieren shows on the teal as well
TEAL_SAND_LIGHT = teal_sand(0.65 / 1.65)


def colormap(name):
    """teal-sand, teal-sand-light or a matplotlib colormap."""
    maps = {"teal-sand": TEAL_SAND, "teal-sand-light": TEAL_SAND_LIGHT}
    return maps[name] if name in maps else matplotlib.colormaps[name]


def schlieren(q, dx, dy=None, lims=(1.0, 100.0), darkest=0.5):
    """Shading in [0, darkest] from |grad q|: none below lims[0], darkest above lims[1], and
    linear in log |grad q| in between. q is indexed [i, j]."""
    gx, gy = np.gradient(q, dx, dy or dx)
    g = np.hypot(gx, gy)
    return darkest * np.clip(np.log(np.maximum(g, 1e-300) / lims[0]) / np.log(lims[1] / lims[0]), 0.0, 1.0)


def _runs(b):
    d = np.diff(np.concatenate([[0], b.astype(np.int8), [0]]))
    return np.where(d == 1)[0], np.where(d == -1)[0]


def level_outlines(fr, lev):
    """Outline of the region covered by the patches of level lev, as line segments in the
    coordinates of the frame. Edges on the boundary of the frame are left out."""
    m = fr.level.T >= lev                                         # [j, i]
    p = np.pad(m, 1, mode="edge")
    x0, _, y0, _ = fr.extent
    X = lambda i: x0 + i * fr.dx
    Y = lambda j: y0 + j * fr.dy
    segs = []
    v = p[1:-1, 1:] != p[1:-1, :-1]                               # vertical edges
    for i in np.where(v.any(0))[0]:
        segs += [[(X(i), Y(a)), (X(i), Y(b))] for a, b in zip(*_runs(v[:, i]))]
    h = p[1:, 1:-1] != p[:-1, 1:-1]                               # horizontal edges
    for j in np.where(h.any(1))[0]:
        segs += [[(X(a), Y(j)), (X(b), Y(j))] for a, b in zip(*_runs(h[j]))]
    return segs


def box_outlines(fr, lev):
    """Edges of every patch box of level lev, cut to the frame, as line segments in the
    coordinates of the frame. Edges on the boundary of the frame are left out, and an edge that
    two boxes share is drawn once."""
    x0, x1, y0, y1 = fr.extent
    tol = 1e-6 * min(fr.dx, fr.dy)
    inside = lambda v, a, b: a + tol < v < b - tol                # noqa: E731
    segs = set()
    for L, a, b, c, d in fr.boxes:
        a, b, c, d = max(a, x0), max(b, y0), min(c, x1), min(d, y1)
        if L != lev or c - a <= tol or d - b <= tol:
            continue
        segs |= {(x, b, x, d) for x in (a, c) if inside(x, x0, x1)}          # vertical edges
        segs |= {(a, y, c, y) for y in (b, d) if inside(y, y0, y1)}          # horizontal edges
    return [[(p, q), (r, t)] for p, q, r, t in sorted(segs)]


def draw_levels(ax, fr, scale=1.0, lw=0.4, alpha=0.5, zorder=1.0, mirror=None, outlines=level_outlines):
    """Outline each AMR level of the frame on ax, thin and light, in the color of the level: the
    region it covers (outlines=level_outlines) or every patch box (box_outlines, as draw_boxes);
    returns the legend handles. mirror = 'x' or 'y' also draws the outlines reflected about the
    lower edge."""
    x0, _, y0, _ = fr.extent
    handles = []
    for L in range(1, int(fr.level.max()) + 1):
        s = np.array(outlines(fr, L)).reshape(-1, 2, 2)
        if mirror == "x":
            s = np.concatenate([s, s * [-1, 1] + [2 * x0, 0]])
        elif mirror == "y":
            s = np.concatenate([s, s * [1, -1] + [0, 2 * y0]])
        ax.add_collection(LineCollection(scale * s, colors=LEVEL_COLORS[L], linewidths=lw, alpha=alpha,
                                         zorder=zorder + 0.1 * L))
        handles.append(Line2D([], [], color=LEVEL_COLORS[L], lw=1.2, alpha=0.85))
    return handles


def draw_boxes(ax, fr, scale=1.0, lw=0.4, alpha=0.5, zorder=1.0, mirror=None):
    """Outline every patch box of each AMR level of the frame on ax in the color of its level;
    the arguments and the legend handles returned are those of draw_levels."""
    return draw_levels(ax, fr, scale, lw, alpha, zorder, mirror, outlines=box_outlines)


# ---------------------------------------------------------------------------
# Layout and output
# ---------------------------------------------------------------------------

def panel_grid(nrows, ncols, aspect, width=10.0, left=0.55, right=0.1, top=0.1, bottom=0.45, wgap=0.12,
               hgap=0.16):
    """Figure with nrows x ncols panels of equal size, height/width = aspect, laid out in inches:
    the margins and gaps are fixed and the height of the figure follows. Returns (fig, axes)."""
    pw = (width - left - right - (ncols - 1) * wgap) / ncols
    ph = pw * aspect
    height = top + nrows * ph + (nrows - 1) * hgap + bottom
    fig = plt.figure(figsize=(width, height))
    axes = np.empty((nrows, ncols), object)
    for r in range(nrows):
        for c in range(ncols):
            axes[r, c] = fig.add_axes([(left + c * (pw + wgap)) / width,
                                       (bottom + (nrows - 1 - r) * (ph + hgap)) / height, pw / width, ph / height])
    return fig, axes


def light_frame(ax):
    """Light gray frame with small outward ticks, for image panels."""
    for s in ax.spines.values():
        s.set_linewidth(0.6)
        s.set_color(FRAME)
    ax.tick_params(labelsize=8, length=2.5, width=0.6, color=FRAME, labelcolor=INK2, pad=2, which="both",
                   direction="out", top=False, right=False)
    ax.minorticks_off()


def colorbar(fig, cax, cmap, norm, label, box=False):
    """Color bar in cax, with the light frame of the image panels or, with box, a dark one."""
    cb = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), cax=cax,
                      ticks=[0, 0.5, 1] if (norm.vmin, norm.vmax) == (0, 1) else None)
    edge = INK if box else FRAME
    cb.set_label(label, fontsize=10 if box else 9, color=INK)
    cb.outline.set_linewidth(0.8 if box else 0.6)
    cb.outline.set_edgecolor(edge)
    cb.ax.tick_params(labelsize=9 if box else 8, length=2.5, width=0.6, color=edge, labelcolor=INK)
    cb.ax.minorticks_off()
    return cb


def save(fig, out, width_px=4800, tight=True, hires=True):
    """Save fig as out, about width_px wide (trimmed to its content when tight), and with
    hires a 400 dpi copy named <out>_400dpi."""
    out = Path(out)
    with plt.rc_context({"savefig.bbox": "tight" if tight else "standard", "savefig.facecolor": "white"}):
        fig.canvas.draw()
        w = fig.get_tightbbox(fig.canvas.get_renderer()).width + 2 * plt.rcParams["savefig.pad_inches"] \
            if tight else fig.get_figwidth()
        fig.savefig(out, dpi=width_px / w)
        if hires:
            fig.savefig(out.with_name(out.stem + "_400dpi" + out.suffix), dpi=400)
    print("wrote", out)


class Movie:
    """H.264 movie (yuv420p) written by ffmpeg from matplotlib figures of the same size, width_px
    wide. The last frame is also saved as <out>_poster.png."""

    def __init__(self, out, fps=12, width_px=1600, crf=20):
        self.out, self.fps, self.width_px, self.crf = Path(out), fps, width_px, crf
        self.proc, self.last, self.n = None, None, 0

    def add(self, fig):
        fig.set_dpi(self.width_px / fig.get_figwidth())
        fig.canvas.draw()
        a = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        a = np.pad(a, ((0, a.shape[0] % 2), (0, a.shape[1] % 2), (0, 0)), constant_values=255)   # even size
        if self.proc is None:
            self.shape = a.shape
            self.proc = subprocess.Popen(
                ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                 "-s", f"{a.shape[1]}x{a.shape[0]}", "-r", str(self.fps), "-i", "-", "-c:v", "libx264",
                 "-pix_fmt", "yuv420p", "-crf", str(self.crf), "-preset", "slow", "-movflags", "+faststart",
                 str(self.out)], stdin=subprocess.PIPE)
        if a.shape != self.shape:
            raise SystemExit("movie frames must all have the same size")
        self.proc.stdin.write(np.ascontiguousarray(a).tobytes())
        self.last, self.n = a, self.n + 1

    def close(self):
        if self.proc is None:
            return
        self.proc.stdin.close()
        if self.proc.wait():
            raise SystemExit("ffmpeg failed")
        poster = self.out.with_name(self.out.stem + "_poster.png")
        plt.imsave(poster, self.last)
        print(f"wrote {self.out} ({self.n} frames, {self.shape[1]}x{self.shape[0]}) and {poster}")


# ---------------------------------------------------------------------------
# Interface measures
# ---------------------------------------------------------------------------

def interface_measures(fr, field="a1"):
    """Extremes of a volume fraction, the length of its 0.5 contour, and its 5 to 95 percent
    thickness, the area where 0.05 < field < 0.95 over that length, in finest cells."""
    import contourpy
    A = fr.q[field]
    lines = contourpy.contour_generator(fr.x, fr.y, A.T).lines(0.5)
    length = sum(np.hypot(*np.diff(l, axis=0).T).sum() for l in lines)
    mixed = ((A > 0.05) & (A < 0.95)).sum() * fr.dx * fr.dy
    return dict(min=A.min(), max=A.max(), length=length, thickness=mixed / length / fr.dx if length else np.nan)


# ---------------------------------------------------------------------------
# Snapshot panels
# ---------------------------------------------------------------------------

def _levels(spec):
    out = []
    for s in spec:
        if ":" in s:
            a, b, n = s.split(":")
            out += list(np.linspace(float(a), float(b), int(n)))
        else:
            out.append(float(s))
    return out


def _fields(a):
    need = [a.field, a.schlieren, a.split, a.contour[0] if a.contour else None]
    return tuple(dict.fromkeys(f for f in need if f))


def resolve_limits(a, frames):
    """Color ranges and schlieren limits not given on the command line, from the frames."""
    def span(k):
        v = np.concatenate([fr.q[k].ravel() for fr in frames])
        return (0.0, 1.0) if k == "a1" else (float(v.min()), float(v.max()))
    if a.field and not a.range:
        a.range = span(a.field)
        if a.field != "a1":
            print(f"--range {a.range[0]:.6g} {a.range[1]:.6g}")
    if a.split and not a.split_range:
        a.split_range = span(a.split)
        print(f"--split-range {a.split_range[0]:.6g} {a.split_range[1]:.6g}")
    if a.schlieren and not a.schlieren_range:
        g1 = max(np.percentile(np.hypot(*np.gradient(fr.q[a.schlieren], fr.dx, fr.dy)), 99.9) for fr in frames)
        a.schlieren_range = (g1 / 300.0, g1)
        print(f"--schlieren-range {g1 / 300.0:.4g} {g1:.4g}")


def _image(fr, field, cmap, norm, a):
    """RGB image [j, i] of field through cmap, white without a field, darkened by the schlieren."""
    rgb = cmap(norm(fr.q[field].T))[..., :3] if field else np.ones(fr.level.T.shape + (3,))
    if a.schlieren:
        rgb = rgb * (1.0 - schlieren(fr.q[a.schlieren], fr.dx, fr.dy, a.schlieren_range).T[..., None])
    return rgb


def draw_panel(ax, fr, a, amr=False):
    """Draw one frame on ax with the drawing options a; returns the AMR legend handles."""
    s = a.scale
    x0, x1, y0, y1 = fr.extent
    kw = dict(origin="lower", interpolation=a.interpolation, zorder=0)
    if a.field and not a.schlieren and not a.mirror:
        ax.imshow(fr.q[a.field].T, cmap=a.cmap, norm=a.norm, extent=np.multiply(fr.extent, s),
                  interpolation_stage="data", **kw)
    elif a.field or a.schlieren:
        ax.imshow(_image(fr, a.field, a.cmap, a.norm, a), extent=np.multiply(fr.extent, s), **kw)
        if a.mirror:
            other = _image(fr, a.split, a.split_cmap, a.split_norm, a) if a.split else \
                _image(fr, a.field, a.cmap, a.norm, a)
            if a.mirror == "y":
                ax.imshow(other[::-1], extent=np.multiply([x0, x1, 2 * y0 - y1, y0], s), **kw)
            else:
                ax.imshow(other[:, ::-1], extent=np.multiply([2 * x0 - x1, x0, y0, y1], s), **kw)
    if a.contour:
        k, lev = a.contour[0], _levels(a.contour[1:])
        ck = dict(levels=lev, colors=INK, linewidths=0.3, zorder=3)
        ax.contour(s * fr.x, s * fr.y, fr.q[k].T, **ck)
        if a.mirror == "y":
            ax.contour(s * fr.x, s * (2 * y0 - fr.y[::-1]), fr.q[k][:, ::-1].T, **ck)
        elif a.mirror == "x":
            ax.contour(s * (2 * x0 - fr.x[::-1]), s * fr.y, fr.q[k][::-1].T, **ck)
    draw = draw_boxes if a.amr_style == "boxes" else draw_levels
    handles = draw(ax, fr, s, lw=a.amr_lw, alpha=a.amr_alpha, zorder=1 if a.frame == "box" else 4,
                   mirror=a.mirror) if amr else []
    xl, yl = a.xlim or (x0, x1), a.ylim or (y0, y1)
    if a.mirror == "x":
        xl = (2 * xl[0] - xl[1], xl[1])
    if a.mirror == "y":
        yl = (2 * yl[0] - yl[1], yl[1])
    ax.set_xlim(s * xl[0], s * xl[1])
    ax.set_ylim(s * yl[0], s * yl[1])
    ax.set_aspect("equal")
    if a.frame == "light":
        light_frame(ax)
    else:                                                  # round tick values: steps of 1, 2 or 5
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(MaxNLocator(nbins="auto", steps=[1, 2, 5, 10]))
    return handles


def snapshot_figure(frames, a):
    """Figure of frames[run][time] with the drawing options a: one row per run and one column per
    time, or with one time and --wrap N the runs in rows of N panels."""
    use_style()
    nruns = len(frames)
    if a.wrap and len(frames[0]) == 1:
        nc = min(a.wrap, nruns)
        nr = -(-nruns // nc)
        cell = lambda r, c: (r * nc + c, 0)                  # noqa: E731  (run, time) of a panel
    else:
        nr, nc = nruns, len(frames[0])
        cell = lambda r, c: (r, c)                           # noqa: E731
    x0, x1, y0, y1 = frames[0][0].extent
    xl, yl = a.xlim or (x0, x1), a.ylim or (y0, y1)
    aspect = (yl[1] - yl[0]) / (xl[1] - xl[0]) * {"y": 2.0, "x": 0.5}.get(a.mirror, 1.0)
    caption_lines = len(a.caption) if a.caption else 0
    top = 0.12 + (0.28 if a.title else 0.0) + 0.19 * caption_lines
    left = (0.62 if a.frame == "box" else 0.55) + (0.3 if a.row_label else 0.0)
    right = 0.95 if a.field else 0.12
    fig, axes = panel_grid(nr, nc, aspect, a.width, left=left, right=right, top=top,
                           bottom=0.52 if a.frame == "box" else 0.45, wgap=0.12, hgap=0.16)
    labels = a.label if a.label is not None else ([] if a.title else ["$t = {t:.3g}$"])
    for r in range(nr):
        for c in range(nc):
            k, j = cell(r, c)
            ax = axes[r, c]
            if k >= nruns:
                ax.set_visible(False)
                continue
            fr = frames[k][j]
            handles = draw_panel(ax, fr, a, amr=a.amr is not None and (not a.amr or k + 1 in a.amr))
            tt = (fr.t - a.t0) * a.tscale
            if a.title and r == 0:
                ax.set_title(a.title.format(t=tt), fontsize=10, color=INK, pad=4)
            txt = labels[min(k, len(labels) - 1)].format(t=tt) if labels else ""
            if txt:
                ax.annotate(txt, (0, 1), xycoords="axes fraction", xytext=(6, -6), textcoords="offset points",
                            ha="left", va="top", fontsize=10, color=INK, zorder=10,
                            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
            if handles:                                   # under the label, in the upper left corner
                below = ScaledTranslation(0, -24 / 72 if txt else 0, fig.dpi_scale_trans)
                leg = ax.legend(handles, [str(L + 1) for L in range(len(handles))], loc="upper left",
                                ncol=len(handles), bbox_to_anchor=(0, 1), bbox_transform=ax.transAxes + below,
                                fontsize=8.5, handlelength=1.3, handletextpad=0.4, columnspacing=0.9,
                                title="AMR level", title_fontsize=8.5, frameon=True, facecolor="white",
                                edgecolor="none", framealpha=0.85)
                leg._legend_box.align = "left"
                leg.set_zorder(10)
            lab = dict(fontsize=9, color=INK, labelpad=1) if a.frame == "light" else {}
            if r == nr - 1 or cell(r + 1, c)[0] >= nruns:
                ax.set_xlabel(a.xlabel, **lab)
                if c < nc - 1 and cell(r, c + 1)[0] < nruns:   # no tick label at the edge next to a panel
                    lo, hi = ax.get_xlim()
                    ax.set_xticks([v for v in ax.get_xticks() if lo <= v < hi - 1e-6 * (hi - lo)])
            else:
                ax.tick_params(labelbottom=False)
            if c == 0:
                ax.set_ylabel(a.ylabel, **lab)
                if r > 0:                                  # no tick label at the edge next to the panel above
                    lo, hi = ax.get_ylim()
                    ax.set_yticks([v for v in ax.get_yticks() if lo <= v < hi - 1e-6 * (hi - lo)])
            else:
                ax.tick_params(labelleft=False)
        if a.row_label and r < len(a.row_label):
            pos = axes[r, 0].get_position()
            fig.text(0.1 / a.width, 0.5 * (pos.y0 + pos.y1), a.row_label[r], rotation=90, ha="left", va="center",
                     fontsize=10, color=INK)
    if a.field:
        box = a.frame == "box"
        t_, b_ = axes[0, -1].get_position().y1, axes[-1, -1].get_position().y0
        xc, wc = axes[0, -1].get_position().x1 + 0.15 / a.width, 0.12 / a.width
        if a.split:
            mid = 0.5 * (t_ + b_)
            g = 0.1 / fig.get_figheight()
            colorbar(fig, fig.add_axes([xc, mid + g, wc, t_ - mid - g]), a.cmap, a.norm, a.colorbar or a.field, box)
            colorbar(fig, fig.add_axes([xc, b_, wc, mid - b_ - g]), a.split_cmap, a.split_norm, a.split, box)
        else:
            colorbar(fig, fig.add_axes([xc, b_, wc, t_ - b_]), a.cmap, a.norm, a.colorbar or a.field, box)
    if a.caption:
        fig.text(axes[0, 0].get_position().x0, 1 - 0.1 / fig.get_figheight(), "\n".join(a.caption), fontsize=8.5,
                 color=INK2, va="top", linespacing=1.5)
    return fig


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _drawing_options(p):
    p.add_argument("--field", default="a1", help="field drawn in color, or none")
    p.add_argument("--cmap", default="teal-sand")
    p.add_argument("--range", type=float, nargs=2)
    p.add_argument("--schlieren", metavar="FIELD")
    p.add_argument("--schlieren-range", type=float, nargs=2, metavar=("G0", "G1"))
    p.add_argument("--contour", nargs="+", metavar="ARG", help="FIELD followed by its levels")
    p.add_argument("--amr", type=int, nargs="*", metavar="RUN")
    p.add_argument("--amr-style", choices=["regions", "boxes"], default="regions",
                   help="with --amr, the region each level covers or every patch box")
    p.add_argument("--amr-lw", type=float, default=0.4, metavar="W", help="width of the AMR lines in points")
    p.add_argument("--amr-alpha", type=float, default=0.5, metavar="A", help="opacity of the AMR lines")
    p.add_argument("--xlim", type=float, nargs=2)
    p.add_argument("--ylim", type=float, nargs=2)
    p.add_argument("--mirror", choices=["x", "y"])
    p.add_argument("--split", metavar="FIELD")
    p.add_argument("--split-cmap", default="viridis")
    p.add_argument("--split-range", type=float, nargs=2)
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--t0", type=float, default=0.0)
    p.add_argument("--tscale", type=float, default=1.0)
    p.add_argument("--title")
    p.add_argument("--label", nargs="+")
    p.add_argument("--row-label", nargs="+")
    p.add_argument("--xlabel", default="$x$")
    p.add_argument("--ylabel", default="$y$")
    p.add_argument("--colorbar")
    p.add_argument("--caption", nargs="+")
    p.add_argument("--wrap", type=int, metavar="N", help="with one time, the runs in rows of N panels")
    p.add_argument("--interpolation", default="antialiased", help="matplotlib interpolation of the images")
    p.add_argument("--frame", choices=["light", "box"], default="light")
    p.add_argument("--width", type=float, default=10.0)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="COMPAS post-processing: plotfiles, snapshot figures and movies",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="See the docstring at the top of this file for the options.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list", help="plotfiles and their times")
    p.add_argument("runs", nargs="+")
    p = sub.add_parser("snapshots", help="figure of panels, one row per run and one column per time")
    p.add_argument("runs", nargs="+")
    p.add_argument("--times", type=float, nargs="+")
    p.add_argument("--out", default="snapshots.png")
    p.add_argument("--width-px", type=int, default=4800)
    _drawing_options(p)
    p = sub.add_parser("movie", help="movie over the plotfiles, one row per run")
    p.add_argument("runs", nargs="+")
    p.add_argument("--tmin", type=float, default=-np.inf)
    p.add_argument("--tmax", type=float, default=np.inf)
    p.add_argument("--every", type=int, default=1, help="use every n-th plotfile")
    p.add_argument("--sync", choices=["time", "span"], default="time",
                   help="runs at the frame time, or each at the same fraction of its own time span")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--out", default="movie.mp4")
    p.add_argument("--width-px", type=int, default=1600)
    _drawing_options(p)
    p = sub.add_parser("interface", help="extremes, contour length and thickness of a volume fraction")
    p.add_argument("runs", nargs="+")
    p.add_argument("--times", type=float, nargs="+")
    p.add_argument("--field", default="a1")
    p.add_argument("--xlim", type=float, nargs=2)
    p.add_argument("--ylim", type=float, nargs=2)
    a = ap.parse_args(argv)
    if a.cmd in ("snapshots", "movie"):
        a.field = None if a.field == "none" else a.field
        if a.split and not a.mirror:
            ap.error("--split needs --mirror")
        if not (a.field or a.schlieren or a.contour):
            ap.error("nothing to draw: give --field, --schlieren or --contour")
        a.cmap, a.split_cmap = colormap(a.cmap), colormap(a.split_cmap)
    return a


def _load(a, path):
    return load(path, _fields(a), a.xlim, a.ylim)


def main(argv=None):
    a = parse_args(argv)
    if a.cmd == "list":
        for run in a.runs:
            print(run)
            for t, p in plotfiles(run):
                print(f"  {t:14.8g}  {p.name}")
    elif a.cmd == "interface":
        w = max(len(run) for run in a.runs)
        print(f"{'run':>{w}s} {'t':>12s} {'min':>11s} {'max - 1':>11s} {'length':>10s} {'thickness':>9s}")
        for run in a.runs:
            for t in (a.times or [None]):
                fr = load(pick(run, t)[1], (a.field,), a.xlim, a.ylim)
                m = interface_measures(fr, a.field)
                print(f"{run:>{w}s} {fr.t:12.6g} {m['min']:11.3e} {m['max'] - 1:11.3e} {m['length']:10.4g} "
                      f"{m['thickness']:9.2f}")
    elif a.cmd == "snapshots":
        frames = []
        for run in a.runs:
            files = [pick(run, t) for t in (a.times or [None])]
            print(run + ": " + ", ".join(f"{p.name} (t = {t:.6g})" for t, p in files))
            frames.append([_load(a, p) for _, p in files])
        resolve_limits(a, [fr for row in frames for fr in row])
        a.norm = Normalize(*a.range) if a.field else None
        a.split_norm = Normalize(*a.split_range) if a.split else None
        save(snapshot_figure(frames, a), a.out, a.width_px)
    elif a.cmd == "movie":
        times = [t for t, _ in plotfiles(a.runs[0], a.tmin, a.tmax)][::a.every]
        if not times:
            raise SystemExit("no plotfiles in the time range")
        files = {run: plotfiles(run, *((a.tmin, a.tmax) if a.sync == "span" else ())) for run in a.runs}
        if not all(files.values()):
            raise SystemExit("a run has no plotfiles in the time range")

        def at(run, i):
            """Plotfile of run for frame i: the one nearest the frame time, or with --sync span
            the one nearest the same fraction of the time span of the run in [tmin, tmax]."""
            f = files[run]
            t = times[i] if a.sync == "time" else f[0][0] + (f[-1][0] - f[0][0]) * i / max(len(times) - 1, 1)
            return min(f, key=lambda g: abs(g[0] - t))[1]

        resolve_limits(a, [_load(a, at(a.runs[0], i)) for i in (0, len(times) - 1)])
        a.norm = Normalize(*a.range) if a.field else None
        a.split_norm = Normalize(*a.split_range) if a.split else None
        movie = Movie(a.out, a.fps, a.width_px)
        for i in range(len(times)):
            fig = snapshot_figure([[_load(a, at(run, i))] for run in a.runs], a)
            movie.add(fig)
            plt.close(fig)
        movie.close()


if __name__ == "__main__":
    main()
