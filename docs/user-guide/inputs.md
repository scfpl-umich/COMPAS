# Run-time options

Run-time options are read from the inputs file, `prob/inputs` in each case, when the run starts,
and changing them needs no rebuild. The options that do are on
[Compile-time options](compile-options.md).

## The inputs file

The file is a list of `name = value` lines, grouped by prefix. `#` starts a comment, and a list
is given as values separated by spaces. Per-direction lists take one value per direction; in 2D a
third value is ignored, so the same line works in 2D and 3D.

```
amr.n_cell      = 64 64 8    # base grid
run.lo_bc       = 0 0 4
FiniteVolume.Scheme = "MUSCL-MC"
```

Any parameter can be overridden on the command line, after the inputs file:

```bash
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs FiniteVolume.Scheme=WENO5 Physics.RiemannSolver=HLL
```

In the tables below, defaults are in bold. An unrecognized value stops the run at startup with a
message that lists the valid values. Parameters marked (AMReX) are read by AMReX, and their
defaults and details are in the [AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/).
The inputs files in `exec/_Tests/` are commented and are a good starting point.

## Run length and domain

| Parameter | Values |
|---|---|
| `max_step` | stop after this many coarse steps (**no limit**) |
| `stop_time` | stop at this physical time (**no limit**) |
| `geometry.prob_lo`, `geometry.prob_hi` | lower and upper corners of the domain, one value per direction (AMReX) |
| `geometry.is_periodic` | `1` periodic, `0` not, one value per direction (AMReX) |
| `geometry.coord_sys` | **`0`**, Cartesian. It is the only value accepted |

The run stops at whichever of `max_step` and `stop_time` is reached first.

## Grid and AMR

| Parameter | Values |
|---|---|
| `amr.n_cell` | base grid, cells per direction |
| `amr.max_level` | number of refinement levels above the base grid (AMReX) |
| `amr.ref_ratio` | refinement ratio between successive levels, one value per level (AMReX) |
| `amr.blocking_factor_x`, `amr.blocking_factor_y`, `amr.blocking_factor_z` | patch sizes are multiples of this (AMReX) |
| `amr.max_grid_size` | largest patch size (AMReX) |
| `amr.n_error_buf` | cells added around tagged cells before patches are built (AMReX) |
| `amr.grid_eff` | fraction of tagged cells a patch must contain (AMReX) |
| `amr.regrid_int` | regrid every N steps (**`2`**) |
| `amr.Interpolater` | coarse–fine interpolation: `pc_interp` (piecewise constant), **`lincc_interp`** (linear), `quartic_interp` |
| `amr.v` | AMReX verbosity (AMReX) |

Where to refine is set under [Refinement](#refinement).

From level 2 up, the ghost cells of a patch are filled from the next coarser level, which covers
only part of the domain. That level must extend $R = \lceil N_g/r\rceil + g$ of its own cells past
the patch. $N_g$ is the number of ghost cells of the scheme: 1 for `GODUNOV`, 2 for the MUSCL
schemes, `WENO3` and `WENO3B`, and 3 for `WENO5` and `WENO5B`. $r$ is the refinement ratio, and $g$
is 0 for `pc_interp`, 1 for `lincc_interp` and 2 for `quartic_interp`. The Phase-Field step always
uses `lincc_interp` for its own ghost cells. AMReX keeps a margin of $a\lceil n_p/a\rceil$ cells,
where $a$ is the blocking factor of the finer level divided by $r$ and $n_p$ is `amr.n_proper`
(AMReX, default 1). Choose the blocking factor so that this margin is at least $R$ in every
direction with a coarse–fine boundary. `exec/_Tests/ShuOsher-5Eq` sets `amr.blocking_factor_x = 8`
for this reason: WENO5 with `lincc_interp` needs a margin of 3, and a blocking factor of 4 gives 2.

## Time stepping

| Parameter | Values |
|---|---|
| `run.cfl` | CFL number (**`0.4`**) |
| `run.TimeIntegrator` | `ForwardEuler`, `TVD-RK2`, **`TVD-RK3`**, `RK4` |
| `run.do_subcycle` | **`1`**: each finer level takes smaller steps, by its refinement ratio. `0`: all levels take the same step |
| `run.do_reflux` | **`1`**: correct the coarse fluxes at coarse–fine boundaries so the update is conservative. `0`: off |
| `run.vnn` | stability number of the viscous and conductive fluxes, $\Delta t =$ `vnn` $\Delta x^2/$(largest diffusivity) (**`0.25`**). Used with `-DDIFFUSION=true`. In the six-equation models the diffusivity is the larger of the momentum and thermal diffusivities of each face ([six-equation fluxes](models.md#viscous-and-conductive-fluxes-in-the-six-equation-models)); in the five-equation models it is $\mu + (\mu_B - \tfrac23\mu) + \kappa$, not divided by $\rho$ or $\rho c_v$ |
| `run.cfl_fast`, `run.cfl_switch` | after `cfl_switch` coarse steps (default `10`) the CFL number becomes `cfl_fast` (default `run.cfl`) |
| `run.timestep_change_limiter` | **`1`**: the time step grows by at most 10% per step. `0`: off |
| `run.check_cfl` | **`1`**: stop if the largest wave speed times the time step exceeds a cell size. `0`: off. CPU builds only |

The time step is shortened to reach `stop_time` exactly, and to land on the output times set by
`amr.t_write_interval` and the output windows ([Output](#output)).

## Boundary conditions

`run.lo_bc` and `run.hi_bc` give one code per direction, for the low and high faces of the
domain.

```
run.lo_bc = 3 3 4      # x, y, z at the low faces
run.hi_bc = 2 2 4      # x, y, z at the high faces
```

| Code | Boundary | Ghost cells |
|---|---|---|
| `0` | interior or periodic | from the periodic image. Used exactly in the directions with `geometry.is_periodic = 1` |
| `1` | inflow | from the case's user boundary function ([Boundary conditions](cases.md#boundary-conditions)), which the case defines |
| `2` | outflow | zero gradient |
| `3` | symmetry | mirror image, with the normal velocity reversed |
| `4` | slip wall | the same as `3` |
| `5` | no-slip wall | mirror image, with all velocity components reversed |

The Phase-Field step uses its own conditions on the volume fractions: periodic in the periodic
directions and zero gradient elsewhere.

## Refinement

A cell is tagged for refinement when any criterion of any listed variable is met.

```
run.refine.ref_vars       = a1 P
run.refine.a1.max_level   = 10
run.refine.a1.value_range = 0.01 0.99  0.01 0.99
run.refine.P.max_level    = 10
run.refine.P.grad         = 0.25 0.25
```

| Parameter | Values |
|---|---|
| `run.refine.ref_vars` | variables used for tagging, any names from [Output variables](#output-variables) |
| `run.refine.<var>.max_level` | cells are tagged on levels below this one (**every level**) |
| `run.refine.<var>.grad` | tag where the jump to a neighboring cell is at least this value |
| `run.refine.<var>.value_greater` | tag where the value is greater than this |
| `run.refine.<var>.value_less` | tag where the value is less than this |
| `run.refine.<var>.value_range` | tag where the value is strictly between a low and a high value |

The thresholds are given per level, one value per level for `grad`, `value_greater` and
`value_less`, and one pair per level for `value_range`. Level 0 comes first. A list must cover
every level that is tagged, from 0 to the lower of `max_level` and `amr.max_level`, minus one, and
a shorter list stops the run at startup. For a region
fixed in space or time, define `spatio_temporal_tag` in the case's `ProblemICBC.H`
([Refinement in space and time](cases.md#refinement-in-space-and-time)).

## Output

| Parameter | Values |
|---|---|
| `amr.case_name` | output folder and prefix. `./plot/MyCase` writes the plotfiles `./plot/MyCase/plt_MyCase00000`, ... and `./plot/MyCase/MyCase.visit` |
| `amr.output_dir`, `amr.plot_file` | alternatives to a path in `amr.case_name`: `amr.output_dir` with `amr.case_name`, or `amr.plot_file` of the form `dir/plt_MyCase`, which also writes to `dir/MyCase/`. The three cannot all be given |
| `amr.plot_int` | a plotfile every N coarse steps (**`-1`**, off). A positive value also writes plotfiles at the start and the end of the run |
| `amr.t_write_interval` | a plotfile every interval of physical time (**`0`**, off) |
| `amr.t_write` | a list of times in increasing order, a plotfile at the first step that reaches each one. Times reached in the same step give one plotfile, and times already passed at a restart are skipped |
| `amr.n_output_windows` | number of time windows with their own output interval (**`0`**) |
| `amr.output_window_<i>.t_lo`, `.t_hi`, `.dt` | window `i`, from `0` to `n_output_windows - 1`: a plotfile every `dt` for `t_lo` $\le t <$ `t_hi`. Outside the windows `amr.t_write_interval` applies |
| `run.output_vars` | variables written to plotfiles ([Output variables](#output-variables)). Without the list, the state variables are written |
| `amr.chk_file` | checkpoint prefix (**`chk`**). `./checkpoints/chk` writes `./checkpoints/chk00500`, ... |
| `amr.chk_int` | a checkpoint every N coarse steps (**`-1`**, off). A positive value also writes one at the start of a new run and one at the end |
| `amr.restart` | checkpoint to restart from, for example `./checkpoints/chk00500` (**empty**, a new run) |
| `run.user_output_int` | call the case's `UserOutputFunction` every N coarse steps (**`0`**, with each plotfile) |
| `run.IO.n_out_files` | number of files each level of a plotfile or checkpoint is written to in parallel (AMReX's default when unset) |
| `run.banner` | **`1`**, `0` to suppress the startup banner |
| `run.timers` | **`1`**: print the wall time of each step and an estimate of the total. `0`: off |
| `run.timestep_sum` | **`1`**: print the sum of state variable `run.sum_var` every step, and stop with a checkpoint if it is NaN. `0`: off |
| `run.sum_var` | index of that state variable, in the order of [Output variables](#output-variables) (**`0`**) |

The cases set `amr.plot_int` to a large number together with `amr.t_write_interval`, which keeps
the plotfiles at the start and the end of the run. `UserOutputFunction` exists only when the case
defines `USER_OUTPUT_FUNC` ([In-situ diagnostics](cases.md#in-situ-diagnostics)).

A run also writes `config_log_<case>.txt` in the output folder, with the build information, the
run-time parameters and copies of the case files.

To restart, give the checkpoint and, if needed, a new stop time:

```bash
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs amr.restart=./checkpoints/chk00500 stop_time=2.0
```

### Output variables

These names are used by `run.output_vars` and `run.refine.ref_vars`. $k$ runs over the phases,
for example `a1rho1 a2rho2 a3rho3`. The state variables come first, in the order of the state
vector, and are written when `run.output_vars` is not given. `rho_w` is a state variable in 3D
only. In 2D, `rho_w` and `w` can still be named, and are zero.

| `PHYSICS` | State variables | Derived variables |
|---|---|---|
| `FIVEEQS` | `a1rho1 a2rho2 rho_u rho_v rho_w rho_E a1` | `rho u v w P T c Mach vorticity_mag q_criterion lambda2` |
| `FIVEEQS_NPHASE` | `a<k>rho<k>`, `rho_u rho_v rho_w rho_E`, `a<k>` | `rho u v w P T c Mach vorticity_mag q_criterion lambda2 Phase PhaseSq` |
| `SIXEQS` | `a1rho1 a2rho2 rho_u rho_v rho_w a1rho1E1 a2rho2E2 a1` | `rho u v w P a1P1 a2P2 T c Mach vorticity_mag q_criterion lambda2` |
| `SIXEQS_IE_NPHASE` | `a<k>`, `a<k>rho<k>`, `a<k>rho<k>e<k>`, `rho_u rho_v rho_w rho_E` | `a<k>p<k> rho u v w P T c Mach vorticity_mag q_criterion lambda2 Phase PhaseSq` |

`P` is the mixture pressure, `T` the temperature, `c` the mixture sound speed and `Mach` the
Mach number. `a1P1` and `a2P2`, and `a<k>p<k>`, are $\alpha_kP_k$. In the N-phase models,
`Phase` is $\sum_k k\,\alpha_k$ and `PhaseSq` is $\sum_k \alpha_k^2$.

`vorticity_mag` is $\lvert\nabla\times\vec u\rvert$, `q_criterion` is
$Q = \tfrac12\left(\lVert\Omega\rVert^2 - \lVert S\rVert^2\right)$ and `lambda2` is the middle
eigenvalue of $S^2 + \Omega^2$, where $S$ and $\Omega$ are the symmetric and antisymmetric parts
of $\nabla\vec u$. They are computed from central differences of $\vec u$, only when requested,
and can also drive refinement through `run.refine.ref_vars`.

## Numerics

| Parameter | Values |
|---|---|
| `FiniteVolume.Scheme` | `GODUNOV`, **`MUSCL`** (minmod), `MUSCL-MC`, `MUSCL-Minmod`, `MUSCL-Superbee`, `MUSCL-vanLeer`, `WENO3`, `WENO3B`, `WENO5`, `WENO5B` |
| `FiniteVolume.Quad` | **`0`**, `1` (quadrature WENO5, `FIVEEQS` in 2D) |
| `FiniteVolume.THINC` | **`0`**, `1` |
| `FiniteVolume.THINC.Beta` | THINC steepness $\beta$ (**`2.3`**) |
| `FiniteVolume.THINC.EPS` | THINC is used only in cells with `EPS` $< \alpha <$ `1 - EPS` (**`1.0e-4`**) |
| `FiniteVolume.THINC.UseNormal` | **`1`**: scale $\beta$ by the interface-normal component in the face direction, and use THINC only on faces where that component exceeds `NormalThreshold`. `0`: use the constant $\beta$ on every face, with no normal check, so `NormalThreshold` is ignored |
| `FiniteVolume.THINC.NormalThreshold` | (**`0.25`**) |
| `FiniteVolume.THINC.Multi` | **`0`**: THINC along each direction. `1`: multi-dimensional THINC: in each THINC cell a $\tanh$ profile across a plane with the interface normal, with the constant $\beta$, is fitted to the cell average, and the face value is its average over the face. `UseNormal` and `NormalThreshold` still select the faces |
| `FiniteVolume.THINC.FaceState` | **`0`**: on a face where THINC is used, both sides take the MUSCL-minmod state for every variable. `1`: both sides keep the state of `FiniteVolume.Scheme`, and only the THINC side is rebuilt from its THINC volume fraction. Use with `PhaseDensity = 1` |
| `FiniteVolume.THINC.PhaseDensity` | **`0`**: the phase densities of a THINC side are the ratio of its face values, $(\alpha_k\rho_k)_f/\alpha_{k,f}$. `1`: they are those of the cell, $\alpha_k\rho_k/\alpha_k$ |
| `FiniteVolume.ID_Bound` | **`0`**, `1` (bound-preserving limiter, `FIVEEQS` and `FIVEEQS_NPHASE` only) |
| `Physics.RiemannSolver` | `LLF`, `HLL`, **`HLLC`** |
| `run.QuadratureIC` | **`1`**: the initial condition is averaged over 4 Gauss–Legendre points per direction in each cell. `0`: the value at the cell center |

`FiniteVolume.Quad = 1` requires `FiniteVolume.Scheme = WENO5` and `FiniteVolume.ID_Bound = 0`,
which the run checks at startup. THINC reconstructs the volume fractions at the interfaces and the
scheme set by `FiniteVolume.Scheme` is used elsewhere. It is available for every model.

`FiniteVolume.THINC.Multi`, `FaceState` and `PhaseDensity` are implemented for `FIVEEQS` and
`FIVEEQS_NPHASE`; the six-equation models stop at startup if any of them is not 0. The defaults,
all three `0`, give the one-dimensional THINC, which works along each direction. For the five-equation
models we recommend all three at `1`. The multi-dimensional face average keeps curved interfaces
from faceting along the grid directions and stops most of the break-up of interfaces nearly
aligned with the grid, and the face state of the base scheme keeps the velocity and pressure from
switching between the base scheme and minmod at THINC faces. Set `PhaseDensity = 1` whenever
`FaceState = 1`: near an interface the ratio of the face values of the base scheme can be far from
the phase density. With `Multi = 1` the default `NormalThreshold = 0.25` keeps curved interfaces
rounder than 0.4. The cases that use these options set them in their `prob/inputs`. Known limits:
sheets and filaments a cell or two wide can still break into beads; with `FaceState = 1` a
strongly sheared interface can carry short waves at the scale of the grid, and a mixed cell where
an interface meets a wall at a sharp corner can keep a spurious low pressure for a short time.

The bound-preserving limiter, `FiniteVolume.ID_Bound = 1`, applies to the five-equation models.
In `SIXEQS` and `SIXEQS_IE_NPHASE` the bounds on $\alpha_k$ and $\alpha_k\rho_k$ come from
`WENO3B` and `WENO5B`.

## Physics

| Parameter | Values |
|---|---|
| `Physics.source_term` | **`0`**, `1` (five-equation models) |
| `Physics.source_term_step_switch`, `Physics.source_term_time_switch` | set `Physics.source_term` to `1` once the coarse step is at least `step_switch` and the time at least `time_switch` (defaults: never, `0.0`). Five-equation models |
| `Physics.pressure_relaxation` | **`0`**, `1` (six-equation models) |
| `Physics.pressure_temperature_relaxation` | **`0`**, `1` (`SIXEQS` only) |

What these terms do is on [Models and equations](models.md).

## Phase-Field

| Parameter | Values |
|---|---|
| `PhaseField.phase_field` | **`0`**, `1` |
| `PhaseField.Mechanism` | **`CAC-Adv`**, `CAC-Adv-NormalFace`, `CAC-VanDerWaals`, `CH-Adv`, `CH-Adv-NormalFace`, `CH-VanDerWaals`, `CDI`, `CDI-NormalFace`, `CDI-CompressionInterpolation`, `ACDI`, `ACDI-NormalFace` |
| `PhaseField.ID_ExplicitRC` | **`0`** implicit, `1` explicit (flux-based mechanisms only) |
| `PhaseField.Eta_Multiplier` | interface width $\eta$ in units of the cell size at `amr.max_level`, the largest over the directions (**`1.0`**) |
| `PhaseField.Eta` | a positive value sets $\eta$ directly and overrides `Eta_Multiplier` (**`0.0`**, unset) |
| `PhaseField.Lambda` | factor on the CAC and Cahn–Hilliard terms (**`0.0`**, which is replaced by `1`) |
| `PhaseField.cfl` | CFL number of the Phase-Field substeps (**`0.4`**) |
| `PhaseField.max_step` | largest number of Phase-Field substeps per time step (**`1000000`**) |
| `PhaseField.Switch_Step`, `PhaseField.Switch_Time` | the Phase-Field step starts once the coarse step and the time reach these values (**`0`**, **`0.0`**) |

After each time step of the model, the Phase-Field terms are advanced over the same interval in
substeps limited by `PhaseField.cfl`. The flux-based mechanisms are `CDI`, `ACDI`, `CH` and their
variants; the `CAC` mechanisms are implicit only (`PhaseField.ID_ExplicitRC = 0`). The mechanisms
and the role of $\eta$ are described under
[Phase-Field mechanisms](models.md#phase-field-mechanisms).

With adaptive mesh refinement (`amr.max_level` above 0), use a `PhaseField.Eta_Multiplier` of
0.25 to 0.5 with the implicit `CAC-Adv` mechanism, as the cases provided that turn Phase-Field on
do. On a single level, values up to the default of 1.0 can be used.

`include/Physics_PhaseField/PhaseField_Parameter.H` also reads `Eta_Factor`,
`ID_DegenerateMobility`, `ID_Potential`, `ID_Weight`, `ID_Flux`, `ID_FluxPhase2Flux` and
`ID_Formulation`, which select variants of the discretization. Leave them at their defaults
unless you are working on the method.

### Linear solver

The implicit Phase-Field treatment solves a linear system with the AMReX multigrid solver
(MLMG). Its parameters have the prefix `LinearSystem`.

| Parameter | Values |
|---|---|
| `LinearSystem.tol_rel`, `LinearSystem.tol_abs` | relative and absolute tolerances (**`1.0e-10`**, **`0.0`**) |
| `LinearSystem.max_iter`, `LinearSystem.max_fmg_iter`, `LinearSystem.fixed_iter` | iteration limits (**`100`**, **`0`**, **`100`**) |
| `LinearSystem.verbose`, `LinearSystem.bottom_verbose` | solver output (**`2`**, **`0`**) |
| `LinearSystem.ID_BottomSolver` | **`0`** BiCGStab, `1` smoother, `2` BiCGStab then CG, `3` CG then BiCGStab, `4` CG |
| `LinearSystem.composite_solve` | **`true`**: solve all levels together. `false`: level by level |
| `LinearSystem.LinOp_maxorder` | largest order of the boundary interpolation of the operator (**`3`**) |
| `LinearSystem.agglomeration`, `LinearSystem.consolidation` | MLMG coarsening options (**`true`**, **`true`**) |
| `LinearSystem.semicoarsening`, `LinearSystem.max_semicoarsening_level` | (**`false`**, **`0`**) |
| `LinearSystem.max_coarsening_level` | (**`30`**) |

The case templates change four of these defaults in `Parm.H`, to `max_iter = 10000`,
`max_fmg_iter = 10000`, `fixed_iter = 10000` and `verbose = 1`. A value in the inputs file
takes precedence over both.

## Equation of state and materials

`prob.EOS.EOS` selects the equation of state for all phases. The forms are given under
[Notation and equations of state](models.md#notation-and-equations-of-state).

| Parameter | Values |
|---|---|
| `prob.EOS.EOS` | `0` ideal gas, `1` stiffened gas, `2` Mie–Grüneisen, `3` Noble-Abel stiffened gas (default set in each case's `Parm.H`) |

The two-phase models, `FIVEEQS` and `SIXEQS`, read one value per phase with the suffixes `_1`
and `_2`:

| Parameter | Symbol | Used by | Default |
|---|---|---|---|
| `prob.EOS.gamma_1`, `prob.EOS.gamma_2` | $\gamma_k$, must be greater than 1 | EOS `0`, `1` and `3` | none, must be given |
| `prob.EOS.pinf_1`, `prob.EOS.pinf_2` | $P_{\infty,k}$ | EOS `1` and `3` | `0.0` |
| `prob.EOS.b_1`, `prob.EOS.b_2` | $b_k$, covolume, must be 0 or positive | EOS `3` | `0.0` |
| `prob.EOS.q_1`, `prob.EOS.q_2` | $q_k$, heat of formation | EOS `3` | `0.0` |
| `prob.EOS.gGamma_1`, `prob.EOS.gGamma_2` | $\Gamma_k$, must be positive | EOS `2` | `0.0` |
| `prob.EOS.pref_1`, `prob.EOS.pref_2` | $p_{\mathrm{ref},k}$ | EOS `2` | `0.0` |
| `prob.EOS.eref_1`, `prob.EOS.eref_2` | $e_{\mathrm{ref},k}$ | EOS `2` | `0.0` |
| `prob.EOS.cv_1`, `prob.EOS.cv_2` | $c_{v,k}$ | temperature | `1.0` |
| `prob.EOS.mu_1`, `prob.EOS.mu_2` | $\mu_k$, shear viscosity | `-DDIFFUSION=true` | `0.0` |
| `prob.EOS.muB_1`, `prob.EOS.muB_2` | $\mu_{B,k}$, bulk viscosity | `-DDIFFUSION=true` | `0.0` |
| `prob.EOS.kappa_1`, `prob.EOS.kappa_2` | $\kappa_k$, thermal conductivity | `-DDIFFUSION=true` | `0.0` |

The N-phase models, `FIVEEQS_NPHASE` and `SIXEQS_IE_NPHASE`, read lists with one value per
phase, for example `prob.EOS.gamma = 1.4 1.6 4.4`. The names are `prob.EOS.gamma`, `pinf`,
`gGamma`, `pref`, `eref`, `b`, `q`, `cv`, `mu`, `muB` and `kappa`. A list that is given must have
exactly `NPHASE` values. A list that is not given is zero, including `cv`, so give `cv` whenever
the temperature is used.

The covolume and the heat of formation (`b_1`, `b_2`, `q_1`, `q_2`, or the lists `b` and `q`)
are read only with EOS `3`. Give them in the units of the case: a case that divides the pressure
scales in its `Parm.H` needs $b$ and $q$ scaled to match. Liquid water, for example, has
$\gamma = 1.19$, $P_\infty = 7.028\times10^{8}$ Pa, $b = 6.61\times10^{-4}$ m³/kg,
$q = -1\,177\,788$ J/kg and $c_v = 3610$ J/(kg K) between 300 and 500 K (Le Métayer and Saurel
2016). With EOS `3`, `FiniteVolume.Quad = 1` is not available: its characteristic decomposition
assumes a stiffened gas.

At startup the run checks that each per-phase list has `NPHASE` values, that each $\gamma_k$ is
greater than 1 with EOS `0`, `1` or `3`, that each $b_k$ is 0 or positive with EOS `3`, and that
each $c_{v,k}$ is positive where the solver needs the temperature: with `-DDIFFUSION=true`, and in
`SIXEQS` with `Physics.pressure_temperature_relaxation = 1`. Its message names the parameter.

Every other `prob.*` parameter belongs to the case and is read by its `Parm.H`
([Parm.H](cases.md#parmh)).

## Convergence studies

With `CONVERGENCE` defined as `true` in the case's `ProblemICBC.H`, COMPAS compares the solution
at the end of the run with the case's exact solution and writes the $L_2$ and
$L_\infty$ errors of every state variable ([Exact solution](cases.md#exact-solution)).

| Parameter | Values |
|---|---|
| `convergence.conv_output_file` | the errors are written to this name with `.0` appended (**`err`**) |
| `convergence.conv_prob_name` | a label written at the top of that file (**empty**) |

`scripts/convergence.py` sets these and runs a case over several resolutions
([Convergence script](testing.md#convergence-script)).
