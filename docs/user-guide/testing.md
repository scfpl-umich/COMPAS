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
| `--only PATTERN` | all cases | Run only the cases whose directory name contains `PATTERN`, and the base runs whose label contains it (only their builds are made) |
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
bash scripts/test_cases.sh --only Advection-6Eq     # base runs whose label contains the pattern
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

1. runs `make realclean`, then builds each build of the case once, before its first base run,
   with `make -jN BUILDARGS ARGS`, N being the number of CPUs unless ARGS sets `-j`; BUILDARGS are
   the make arguments of the base run (`BUILD=6Eq`, see
   [Several builds and base runs in one case](#several-builds-and-base-runs-in-one-case)), none for
   a case with one build,
2. runs the executable of each base run with its inputs file (`prob/inputs` for a case with one
   build),
3. runs it again for each variant line of that base run in `variants.txt`,
4. removes the output and every build.

Before building, the runner asks make for the executable of each build
(`make BUILDARGS print-executable`) and refuses two builds that would write the same executable. In
a case with several builds it deletes `tmp_build_dir` after each build, so the object files of only
one build are on disk at a time; the executables stay until the end of the case.

Before each run, except the `+` variants described below, it calls `exec/clean.sh` in the case
directory, which deletes and recreates `plot/` and `checkpoints/`.

:::{warning}
Running the suite deletes the plotfiles and checkpoints in the case directories of
`exec/_Tests/`. Copy any output you want to keep before running it.
:::

## Variants

`variants.txt` lets one executable cover many run-time options. Each line is a variant name
followed by overrides appended to the command line of the base run, and text after `#` is a
comment. From the `[Advection-5Eq]` section of `exec/_Tests/Advection/variants.txt`:

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

## Several builds and base runs in one case

A case directory can hold several builds, one per model, and several base runs. `exec/_Tests/Sod`,
for example, builds the Sod shock tube with each of the four models. Its `GNUmakefile` maps `BUILD`
to the DEFINES of each build, and its `prob/` serves every build, with `#if (PHYSICS == ...)`
where the models differ:

| Command | Builds | Executable |
|---|---|---|
| `make` | the default, the first of `BUILDS` in the `GNUmakefile` (here `BUILD=5Eq`) | `main3d.gnu.MPI.ex` |
| `make BUILD=6Eq` | another build | `main3d.gnu.MPI.6Eq.ex` |

The default build keeps the plain executable name and its first base run the plain `prob/inputs`, so
`make && ./main2d.gnu.MPI.ex prob/inputs` (`main3d` in a 3D case) works in every case directory.
The other builds add their name to the executable and the object directory
(`USERSuffix := .$(BUILD)`), so the builds can sit side by side.

In `variants.txt` a line `[LABEL] ARGS` starts a base run named `LABEL`: the words of ARGS are
passed to `make` (`BUILD=...`), except `inputs=FILE`, its inputs file. The variant lines below it
belong to it, and base runs with the same `BUILD` share one build. The labels are the names the
runs had as separate cases. From `exec/_Tests/Diffusion/variants.txt`:

```
[Couette2Layer-6Eq]          BUILD=6Eq     inputs=prob/inputs.Couette2Layer-6Eq
mid-cell                max_step=20 prob.h=0.515625
no-relaxation           max_step=20 Physics.pressure_relaxation=0

[ViscousShockTube-6Eq]       BUILD=6Eq     inputs=prob/inputs.ViscousShockTube-6Eq
relax-pT                max_step=10 Physics.pressure_relaxation=0 Physics.pressure_temperature_relaxation=1
```

The two base runs share the `6Eq` build; `prob.problem` in each inputs file chooses the problem
(the list is in `exec/_Tests/Diffusion/prob/Parm.H`). A `variants.txt` without such lines
describes one build, run with `prob/inputs` and named after its directory, as in
`RichtmyerMeshkov-5Eq`. Labels must be unique in the suite; the logs and the summary use them.

To build and run one base run by hand, take its `BUILD` and inputs file from `variants.txt`:

```bash
cd exec/_Tests/Diffusion
grep '^\[' variants.txt                       # the base runs: label, BUILD and inputs file
make -j4 BUILD=6Eq                            # main2d.gnu.MPI.6Eq.ex
mpirun -n 4 ./main2d.gnu.MPI.6Eq.ex prob/inputs.Couette2Layer-6Eq
make BUILD=6Eq cleanconfig                    # removes this build only
make realclean                                # removes every build of the directory
```

`make BUILD=6Eq cleanconfig` deletes the executable and the object files of that build and leaves
the others; `make realclean` (and `make clean`) deletes `tmp_build_dir` and every executable of the
directory. The runs of one directory share `plot/` and `checkpoints/`, and `exec/clean.sh` empties
both, so run one base run at a time per directory. The configuration log of a run copies the
inputs file the run reads (`prob/inputs` or `prob/inputs.<label>`).

## What the cases cover

The 14 case directories hold 40 base runs in 31 builds; with their variants they make 275 runs.

| Directory | Builds | Base runs (labels) | What it tests |
|---|---|---|---|
| `Advection` | 5Eq, 5Eq-N, 6Eq, 6Eq-N | `Advection-5Eq`, `Advection-5Eq-N`, `Advection-6Eq`, `Advection-6Eq-N`, `Advection-THINC-5Eq`, `Advection-PF-5Eq` | Four shapes advected at uniform pressure and velocity. Each model runs every reconstruction scheme, Riemann solver and equation of state, and each family of Phase-Field mechanisms; `Advection-5Eq` also runs all eleven mechanisms and the time integrators, and the `eos-nasg-b0` variants check that the NASG EOS reproduces the stiffened gas with $b = q = 0$. `Advection-THINC-5Eq` and `Advection-PF-5Eq` are inputs files of the 5Eq build with THINC and with Phase-Field on, to try each out of the box |
| `Diffusion` | 5Eq, 5Eq-N, 6Eq, 6Eq-N, 6Eq-N3 | `Couette2Layer-*`, `HeatConduction2Mat-*`, `ViscousShockTube-*` (5Eq, 6Eq, 6Eq-N), `Advection-Viscous-6Eq`, `-6Eq-N`, `-5Eq-N` | The viscous and conductive fluxes (`-DDIFFUSION=true`): a two-layer Couette flow with a moving wall and conduction between two materials, across a material interface; the viscous shock tube; viscous advection of a drop (with an absent third phase in the N-phase model, build 6Eq-N3) and of five materials |
| `Sod` | 5Eq, 5Eq-N, 6Eq, 6Eq-N | `Sod-5Eq`, `Sod-5Eq-N`, `Sod-6Eq`, `Sod-6Eq-N` | The 3D Sod shock tube with AMR along x, y and z, the boundary conditions, and the two-material water-air shock tube |
| `NASG` | 5Eq, 5Eq-N, 6Eq, 6Eq-N | `NASG-5Eq`, `NASG-5Eq-N`, `NASG-6Eq`, `NASG-6Eq-N` | The Noble-Abel stiffened gas: the air-water shock tube, a water shock tube and a water drop |
| `ShuOsher-5Eq`, `ShockVortex-5Eq`, `IsentropicVortex-5Eq` | one each | the directory name | The reconstruction schemes on the Shu-Osher tube, a vortex through a stationary shock with AMR, and the convergence case with WENO5 |
| `StaticDrop` | 5Eq, 5Eq-N, 6Eq, 6Eq-N | `StaticDrop-5Eq`, `StaticDrop-5Eq-N`, `StaticDrop-6Eq`, `StaticDrop-6Eq-N` | Surface tension (`-DSURFACE_TENSION=true`) on each model, with the initial pressure projection, both energy forms, Phase-Field, THINC, the stress without smoothing and AMR |
| `COMPAS-STL-5Eq`, `RichtmyerMeshkov-5Eq`, `RichtmyerMeshkov-Multimode-5Eq`, `Jet-Inflow-5Eq`, `NonsphericalCollapse-6Eq` | one each | the directory name | One user hook or feature each: the STL reader, the user output with checkpoint and restart, run-time static GPU arrays, a user boundary condition, user refinement tagging in 3D |
| `RayleighTaylor` | 2D, 3D (`BUILD` sets `DIM`) | `RayleighTaylor-5Eq`, `RayleighTaylor-3D-5Eq` | A user source term, in 2D and 3D |

## Pass and fail

Each run gets one status, checked in this order:

| Status | Meaning | Result |
|---|---|---|
| `NAN` | The run log contains `nan` as a word, in any case | fail |
| `TIMEOUT` | The run was still going at the time limit | pass |
| `RUN_FAIL` | The run exited with a nonzero code | fail |
| `UNUSED_PARAM` | AMReX listed an override among its unused ParmParse variables, so the variant tested nothing | fail |
| `OK` | The run finished | pass |

A base run whose build fails is recorded as `FAILED` and its runs are skipped, and a build that
produces no executable, or the executable of another build of the case, is recorded as `NOEXE`.
Both count as failures.

AMReX lists the unused variables when the run ends, so the `UNUSED_PARAM` check needs a run that
finishes before the time limit. This is why the variants set a small `max_step`.

The runner exits with 0 if every run passed, 1 if any failed, and 2, before building anything, if
`--only` matched no case, an option is not recognized, a label appears twice in the suite, or a
`variants.txt` has a variant line before its first `[LABEL]` line. The build logs are `scripts/logs/<label>.log`, named after the
first base run of the build, the run logs `scripts/logs/<label>_run.log` and
`scripts/logs/<label>__<variant>_run.log`, and the summary is `scripts/test_summary.csv`, with one
row per run, the first column being `<label>` or `<label>:<variant>`:

```
case,build_status,run_status,build_seconds,run_seconds,total_seconds
```

## Adding a test

A new case in `exec/_Tests/` is picked up automatically. Add its name to `test_list.txt` to fix
its place in the order and remove the warning. A new run-time option is tested by a line in the
`variants.txt` of a case for each model it applies to. A new model of an existing problem is a new
`BUILD` in the case's `GNUmakefile`, `#if (PHYSICS == ...)` branches in `prob/` where the model
differs, an inputs file and a `[LABEL] BUILD=... inputs=...` section. A new problem of the
`Diffusion` case is a header `prob/<Problem>.H` with its functions in a namespace of its name, its
parameters and a value of `DiffusionProblem` in `Parm.H`, a line in each dispatch of
`ProblemICBC.H`, and a section per model. What a pull request needs is on the
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
