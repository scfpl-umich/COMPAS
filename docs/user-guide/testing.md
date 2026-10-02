# Testing

The test suite builds every case in `exec/_Tests/` and runs it briefly, once with its inputs file
and once for each of its variants, to check that each model and each option builds and runs
stably. The comparisons with reference solutions are on the
[Verification and validation](../verification.md) page.

## Running the suite

From the repository root,

```bash
bash scripts/test_cases.sh
```

`scripts/test_cases.py` is an equivalent runner that needs only Python 3, with the same options:

```bash
python3 scripts/test_cases.py
```

| Option | Default | Meaning |
|---|---|---|
| `--only PATTERN` | all cases | Run only the cases whose directory name contains `PATTERN` |
| `--no-variants` | off | Run only the base inputs of each case |
| `--ranks=N` | 1 | MPI ranks per run. With 1, the executable runs without `mpirun` |
| `--timeout=S` | 10 | Time limit of each run, in seconds |
| `-q`, `--quiet` | off | Print only the progress and the summary. Without it the runners also print the commands, log files and timings of each step |
| `--tail=N` | 40 | Lines of the log printed when a build or run fails |
| `--dont-clean` | off | Keep the output of the last run of each case |
| `-h`, `--help` | | Print the options |
| `-- ARGS` | | Pass everything after `--` to `make`. The build uses one job per CPU unless ARGS sets `-jN` |

The options that take a value accept both `--opt=N` and `--opt N`. An unknown option stops the
runner with an error.

For example,

```bash
bash scripts/test_cases.sh --only Advection-6Eq     # cases whose name contains the pattern
bash scripts/test_cases.sh --no-variants            # base inputs only
bash scripts/test_cases.sh --ranks=4 --timeout=30
bash scripts/test_cases.sh --quiet -- -j4           # arguments after -- go to make
```

The shell runner enforces the time limit with `timeout`, or `gtimeout` from GNU coreutils on
macOS (`brew install coreutils`); without either, it warns once and runs without a limit. The
Python runner has its own timer.

## What the suite does

The cases run in the order of `exec/_Tests/test_list.txt`. A case directory with a `GNUmakefile`
that is not in the list runs after the listed ones, with a warning.

For each case the runner

1. runs `make realclean` and builds the case with `make -jN ARGS`, N being the number of CPUs unless ARGS sets `-j`, using its own `GNUmakefile`,
2. runs the executable with `prob/inputs`,
3. runs it again for each line of the case's `variants.txt`, if there is one,
4. removes the output and the build.

Before each run, except the `+` variants described below, it calls `exec/clean.sh` in the case
directory, which deletes and recreates `plot/` and `checkpoints/`.

:::{warning}
Running the suite deletes the plotfiles and checkpoints in the case directories of
`exec/_Tests/`. Copy any output you want to keep before running it.
:::

## Variants

`variants.txt` lets one executable cover many run-time options. Each line is a variant name
followed by overrides appended to the command line of the base run, and text after `#` is a
comment. From `exec/_Tests/Advection-5Eq/variants.txt`:

```
scheme-weno5            max_step=10 FiniteVolume.Scheme=WENO5
riemann-hll             max_step=10 Physics.RiemannSolver=HLL
pf-acdi                 max_step=10 PhaseField.phase_field=1 PhaseField.Mechanism=ACDI PhaseField.ID_ExplicitRC=1
```

The output is cleaned before every variant, except that a variant whose name starts with `+`
continues from the output of the previous line. `exec/_Tests/RichtmyerMeshkov-5Eq` uses this to
test a restart:

```
restart-write           max_step=5 amr.chk_int=5
+restart-read           max_step=10 amr.restart=./checkpoints/chk00005
```

`Advection-5Eq`, `Advection-5Eq-N`, `Advection-6Eq` and `Advection-6Eq-N` run every
reconstruction scheme, Riemann solver and equation of state, and each family of Phase-Field
mechanisms, on their model. `Advection-5Eq` also runs all eleven mechanisms and the time
integrators, and the four `Sod-*` cases run the boundary conditions and the two-material water-air
shock tube on each model. The `ShearDecay-*`, `HeatConduction-*`, `Advection-Viscous-*` and
`ViscousShockTube-*` cases run the viscous and conductive fluxes (`-DDIFFUSION=true`). The four
`NASG-*` cases run the Noble-Abel stiffened gas on each model, and the `eos-nasg-b0` variants of the
`Advection-*` cases check that it reproduces the stiffened gas with $b = q = 0$. The 34 cases and
their variants make 243 runs.

## Pass and fail

Each run gets one status, checked in this order:

| Status | Meaning | Result |
|---|---|---|
| `NAN` | The run log contains `nan` as a word, in any case | fail |
| `TIMEOUT` | The run was still going at the time limit | pass |
| `RUN_FAIL` | The run exited with a nonzero code | fail |
| `UNUSED_PARAM` | AMReX listed an override among its unused ParmParse variables, so the variant tested nothing | fail |
| `OK` | The run finished | pass |

A case whose build fails is recorded as `FAILED` and its runs are skipped, and a build that
produces no `*.ex` file is recorded as `NOEXE`. Both count as failures.

AMReX lists the unused variables when the run ends, so the `UNUSED_PARAM` check needs a run that
finishes before the time limit. This is why the variants set a small `max_step`.

The runner exits with 0 if every run passed, 1 if any failed, and 2 if `--only` matched no case
or an option is not recognized. The build logs are `scripts/logs/<case>.log`, the run logs `scripts/logs/<case>_run.log` and
`scripts/logs/<case>__<variant>_run.log`, and the summary is `scripts/test_summary.csv`, with one
row per run:

```
case,build_status,run_status,build_seconds,run_seconds,total_seconds
```

## Adding a test

A new case in `exec/_Tests/` is picked up automatically. Add its name to `test_list.txt` to fix
its place in the order and remove the warning. A new run-time option is tested by a line in the
`variants.txt` of a case for each model it applies to. What a pull request needs is on the
[Contributing](../about/contributing.md) page.

## Convergence script

`scripts/convergence.py` measures the order of accuracy of a case that defines `CONVERGENCE` and
`exact_solution` in its `ProblemICBC.H` (see [Setting up a case](cases.md#exact-solution)). With
them, COMPAS writes the $L_2$ and $L_\infty$ norms of the error of each conservative variable
on level 0 at the end of the run, to the file set by `convergence.conv_output_file`.

Run the script from the case directory in `exec/_Tests/`:

| Option | Default | Meaning |
|---|---|---|
| `--run` | | Run the case at each resolution and save the errors |
| `--plot` | | Plot the saved errors on log-log axes |
| `--res N ...` | `8 16 32 64 128 256` | Resolutions; each run sets `amr.n_cell=N N 1` |
| `--norm NORM` | `Linf` | `Linf` or `L2`, the norm saved by `--run` and named on the plot |
| `--nprocs N` | 1 | MPI ranks, through `mpirun -n N`. With 1, the executable runs without `mpirun` |
| `--exe EXE` | the only `*.ex` in the case directory | Executable |
| `--inputs FILE` | `./prob/inputs` | Inputs file |
| `--vars VAR ...` | all | Variables to plot. `--run` always saves every variable |
| `--csv FILE` | `convergence_all.csv` | File written by `--run` and read by `--plot` |
| `--clear true\|false` | `true` | Run `exec/clean.sh` in the case directory before the runs |

`--run` saves the chosen norm of the error of every variable against the grid spacing
$h = (x_\mathrm{hi} - x_\mathrm{lo})/N$, with the domain read from `geometry.prob_lo` and
`geometry.prob_hi` in the inputs file, and prints the order between each pair of successive
resolutions. It stops with an error if a run writes no error file or a file for another
resolution. `--plot` draws each variable with a fitted line through the last two points.

Outside a virtual environment, the script sets up and uses the one described under
[Python environment](../getting-started/install.md#python-environment).

`exec/_Tests/IsentropicVortex-5Eq` is the convergence case. The commands are in
[Order of accuracy](../verification.md#order-of-accuracy).

## Post-processing script

`scripts/postprocess.py` reads the plotfiles of a 2D run with yt and draws them in the style of
the figures on these pages. A run is given by its plot directory, for example `plot/KH-M0.4`, and
the fields are composited on the grid of the finest level, each level interpolated linearly and
used where no finer level covers it. Run it from the case directory:

| Command | What it does |
|---|---|
| `list RUN ...` | Print the plotfiles of each run and their times |
| `snapshots RUN ... [--times T ...]` | One row of panels per run and one column per time, at the nearest plotfiles, or at the last one without `--times`. Writes a figure about 4800 px wide and a 400 dpi copy |
| `movie RUN ... [--tmin T] [--tmax T]` | The same panels over every plotfile of the first run in the time range, one row per run, as an H.264 MP4, 1600 px wide by default, with the last frame as a poster image. With `--sync span`, runs that end at different times all run from start to end together. Needs `ffmpeg` |
| `interface RUN ... [--times T ...]` | The extremes of a volume fraction, the length of its 0.5 contour, and its 5 to 95 percent thickness in finest cells |

The main drawing options of `snapshots` and `movie` are

| Option | Meaning |
|---|---|
| `--field NAME` | Field drawn in color, `a1` by default, or `none` |
| `--cmap NAME` | `teal-sand`, the volume-fraction map of these pages, its lighter part `teal-sand-light`, or a matplotlib colormap |
| `--schlieren FIELD` | Darken the colors with a light schlieren of the gradient of `FIELD` |
| `--contour FIELD LEVELS` | Contour lines; `START:STOP:N` gives `N` equally spaced levels |
| `--amr [N ...]` | Outline the region covered by each AMR level, in every panel or only in the panels of runs number `N`, counted from 1 |
| `--amr-style regions\|boxes` | With `--amr`, the outline of the region each level covers (the default), or every patch box of each level, in the color of its level |
| `--wrap N` | With one time, lay the runs out in rows of `N` panels instead of one row per run |
| `--xlim X0 X1`, `--ylim Y0 Y1` | The window |
| `--mirror x\|y`, `--split FIELD` | Reflect each panel about its lower edge, a symmetry plane, and draw another field in the reflected half |
| `--title FMT`, `--label FMT` | Column titles and panel labels, with the time as `{t}` |

The docstring at the top of the script lists every option. For example, from
`exec/KelvinHelmholtz-5Eq`,

```bash
python3 ../../scripts/postprocess.py list plot/KH-M0.4
python3 ../../scripts/postprocess.py snapshots plot/KH-M0.4 --times 65 75 85 95 --ylim -8.75 8.75
python3 ../../scripts/postprocess.py movie plot/KH-M0.4 --ylim -8.75 8.75 --fps 8 --out kh.mp4
```

The analysis scripts of the Kelvin-Helmholtz, shock refraction, double Mach reflection and water
hammer cases import the same functions for the figures that need more than these options, and
hold only what belongs to their case, such as an exact solution or a measurement. Outside a
virtual environment, the script and these analysis scripts set up and use the environment
described under [Python environment](../getting-started/install.md#python-environment).

The double Mach video on the
[Verification and validation](../verification.md#double-mach-reflection) page is drawn with `movie`
and `--width-px 2400`, from a run that writes plotfiles more often.
