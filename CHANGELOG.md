# Changelog

Changes that users need to know about, newest first. A change that forces an edit to a case goes
under "Case interface", with a short before-and-after example.

## Unreleased

First public release.

### Added

- Viscous and conductive fluxes for the six-equation models: `-DDIFFUSION=true` now builds
  with `PHYSICS=SIXEQS` and `PHYSICS=SIXEQS_IE_NPHASE`. Each phase has its own stress and its
  own heat flux with its own temperature (Huang, arXiv:2609.18085), and mass, momentum and total
  energy stay conservative, also across AMR levels. The face viscosity is the harmonic mean of
  the mixture viscosities of the two cells, and with the pressure-temperature relaxation of
  `SIXEQS` the face conductivity is the harmonic mean of the cell conductivities. Without that
  relaxation, and always in `SIXEQS_IE_NPHASE`, heat does not pass from one phase to another (such
  a run prints a warning when a conductivity is nonzero). See
  [Viscous and conductive fluxes in the six-equation models](docs/user-guide/models.md#viscous-and-conductive-fluxes-in-the-six-equation-models).
  The previous `SIXEQS_IE_NPHASE` diffusion kernels, which did not compile with the core, are
  replaced. With `DIFFUSION=true`, `SIXEQS` checks that `prob.EOS.cv_1` and `prob.EOS.cv_2` are
  positive.
- Noble-Abel stiffened-gas equation of state (Le Métayer and Saurel 2016), `prob.EOS.EOS = 3`,
  for all four models. The covolume and the heat of formation are new keys, `prob.EOS.b_1`,
  `b_2`, `q_1`, `q_2` for the two-phase models and the lists `prob.EOS.b`, `prob.EOS.q` for the
  N-phase models (default 0), so no case file changes. With `b = q = 0` it reproduces the
  stiffened gas. A state with $1 - b_k\rho_k \le 0$ stops the run. New test cases `NASG-5Eq`,
  `NASG-5Eq-N`, `NASG-6Eq`, `NASG-6Eq-N`, and `eos-nasg-b0` variants of the `Advection-*` cases.
- Surface tension for all four models, built with `-DSURFACE_TENSION=true` (default `false`,
  which leaves the code unchanged). The capillary stress
  $\sum_k\sigma_k(\lvert\nabla\alpha_k\rvert\mathbb{I} - \nabla\alpha_k\otimes\nabla\alpha_k/\lvert\nabla\alpha_k\rvert)$
  is a conservative momentum face flux, refluxed across AMR levels, computed from the volume
  fractions smoothed by a 1-2-1 filter (`SurfaceTension.smoothing`). The energy gets either the
  work of the capillary force (`SurfaceTension.energy_form = work`, the default) or a conservative
  energy flux (`conservative`). The capillary time step limit is applied, and
  `SurfaceTension.init_pressure = 1` projects the initial pressure onto the discrete capillary
  force. New keys `SurfaceTension.sigma` (pairwise coefficients for the N-phase models),
  `energy_form`, `init_pressure`, `init_pressure_weight` and `smoothing`, read by the models'
  physics parameters, so no case file changes. See [Surface tension](docs/user-guide/models.md#surface-tension).
  New test cases `StaticDrop-5Eq`, `StaticDrop-5Eq-N`, `StaticDrop-6Eq` and `StaticDrop-6Eq-N`.
  Experimental: the spurious currents grow as the grid is refined, and with AMR the whole
  interface must stay on one level (see the user guide).
- Test cases (labels of the `Diffusion` directory, see "Changed") `ViscousShockTube-6Eq`,
  `ViscousShockTube-6Eq-N`, `Advection-Viscous-6Eq` and `Advection-Viscous-6Eq-N`, and the
  material-interface tests
  `Couette2Layer-5Eq`, `Couette2Layer-6Eq`, `Couette2Layer-6Eq-N` (two-layer Couette flow, with a
  moving wall through the user boundary condition) and `HeatConduction2Mat-5Eq`,
  `HeatConduction2Mat-6Eq`, `HeatConduction2Mat-6Eq-N` (conduction between two materials), with
  a `mid-cell` variant of `Couette2Layer-6Eq` (the interface inside a cell). `Sod-6Eq` has a
  `water-air-512` variant (the water-air shock tube on 512 uniform cells).
- `Physics.pressure_relaxation_trace` (six-equation models, default `0`): a phase with
  $0 < \alpha_k \le$ this value keeps its volume fraction in the pressure relaxation and reaches
  the common pressure through its energy. See
  [Safeguarded pressure relaxation](docs/user-guide/models.md#safeguarded-pressure-relaxation).

### Changed

- Test suite (`exec/_Tests`): one directory per problem family. `Advection`, `Diffusion`, `Sod`,
  `NASG`, `StaticDrop` and `RayleighTaylor` hold several builds each, one per model
  (one per dimension in `RayleighTaylor`), chosen by `BUILD` in the directory's `GNUmakefile`, and
  one `prob/` with `#if (PHYSICS == ...)` where the models differ; in `Diffusion`, `prob.problem`
  in the inputs file picks the problem. The default build of a directory keeps the plain
  executable name and `prob/inputs`; the others are `make BUILD=6Eq` and so on, with their own
  executable (`main2d.gnu.MPI.6Eq.ex`). In `variants.txt` a line `[LABEL] BUILD=... inputs=...`
  starts each base run, and the labels are the former case names, so `--only`, the logs and the
  summary keep them. The release cases moved: `Advection-5Eq`, `Advection-5Eq-N`, `Advection-6Eq`,
  `Advection-6Eq-N`, `Advection-THINC-5Eq` and `Advection-PF-5Eq` to `Advection`; `Sod-*` to `Sod`;
  `RayleighTaylor-5Eq` and `RayleighTaylor-3D-5Eq` to `RayleighTaylor`; `ViscousShockTube-5Eq` and
  `Advection-Viscous-5Eq-N` to `Diffusion`. `Advection-THINC-5Eq` and `Advection-PF-5Eq` are now
  inputs files of the 5Eq build of `Advection`: they start from the initial condition of
  `Advection-5Eq` and write `advection.csv` instead of `interface.csv`, with the interface
  thickness `N_I` as its last column, so the three methods can be compared.
  `ViscousShockTube-5Eq` and `IsentropicVortex-5Eq` now run WENO5 (`FiniteVolume.Quad = 0`), as
  the six-equation viscous shock tubes do; the plot folder of `IsentropicVortex-5Eq` is
  `plot/IsentropicVortex`.
  `scripts/test_cases.sh` and `scripts/test_cases.py` build each build of a directory once and
  refuse two builds that would write the same executable; a directory without sections runs as
  before.
  The suite has 14 directories, 31 builds, 40 base runs and 275 runs.
- AMR without subcycling (`run.do_subcycle = 0`): the levels now take every Runge-Kutta stage
  together, and with `run.do_reflux = 1` the coarse faces covered by a finer level take, at every
  stage and before any level is updated, the area-weighted average of its face fluxes: the
  conservative fluxes with their viscous and capillary parts, the face quantities of the
  non-conservative terms, and in the Phase-Field step its fluxes. The update is conservative at
  every stage, with no reflux after the step, and the volume fraction stays within $[0, 1]$ at a
  sharp interface that crosses a coarse-fine boundary, also with `TVD-RK2`, `TVD-RK3` and `RK4`.
  With `run.do_reflux = 0` the stages run the same way without the averaged faces. Subcycled runs
  and runs on one level are unchanged.

### Fixed

- The configuration log `config_log_<case>.txt` copied `./prob/inputs` whatever inputs file the
  run was given. It now copies the inputs file named on the command line. Results are unchanged.
- OpenMP (`USE_OMP = TRUE`): threads could scale a shared face flux twice, or read it after
  another thread had scaled it, and two threads could update the largest wave speed of a cell at
  once, so results depended on the thread timing. Several threads now give the same results as
  one, and builds without OpenMP are unchanged. Runs with an AMReX linear solve (the implicit
  Phase-Field step, `PhaseField.ID_ExplicitRC = 0`,
  and `SurfaceTension.init_pressure = 1`)
  still depend on the number of threads through the order of AMReX's sums, unless
  `amrex.regtest_reduction = 1`.
- `IsentropicVortex-5Eq` and other cases with `CONVERGENCE` did not compile with
  `USE_OMP = TRUE`. Results are unchanged.
- Phase-Field step without subcycling (`PhaseField.phase_field = 1`, `run.do_subcycle = 0`): each
  Phase-Field step also advanced the step counter of the flow, which doubled the plotfile numbers
  and the counts of `amr.regrid_int` and `amr.chk_int`. The Phase-Field steps no longer count, as
  with subcycling.
- Six-equation models with a relaxation and AMR (with `run.do_reflux = 1`): the reflux corrected
  the coarse cells next to a finer level after the relaxation, so those cells ended every step out
  of equilibrium. The reflux now relaxes again the cells it changes. Runs without relaxation,
  without AMR levels or without reflux are unchanged.
- `SIXEQS`, HLLC Riemann solver: the star state now scales the volume fraction like the partial
  densities, so the first-order update of the volume fraction is upwind and keeps it bounded.
- Pressure relaxation of the six-equation models (`Physics.pressure_relaxation = 1`, and
  `Physics.pressure_temperature_relaxation = 1` in `SIXEQS`): a phase energy below the minimum of
  its equation of state, a volume fraction just outside $[0, 1]$, or a root next to an asymptote
  of the relaxation function made it fail (in `SIXEQS_IE_NPHASE` with NaN, which stopped the
  water-air shock tube of `Sod-6Eq-N` at the first step). Where their own solution fails, both
  models now use a safeguarded relaxation that keeps the masses, the momentum and the total energy;
  elsewhere results are unchanged. No relaxation clips or rescales a volume fraction.
  `Sod-6Eq-N@water-air` uses the water of the other `Sod` cases again.
  `Physics.pressure_relaxation_implicit`, read but never passed to the relaxation, is removed.
- `SIXEQS`, `Physics.pressure_temperature_relaxation = 1`: the root of the quadratic of the
  closed form lost its digits for a pressure small against $P_{\infty,k}$. It is now computed in
  a stable form, and with the Noble-Abel stiffened gas the closed form is used only where every
  phase with mass stays below its covolume limit (otherwise the pressure relaxation). Runs with
  this relaxation change by round-off, and more where the closed form had lost its digits.
- `SIXEQS_IE_NPHASE`: the Phase-Field flux kernel left its non-conservative face components
  uninitialized. They were never used, so results are unchanged; the kernel now sets them to 0,
  as the other models do.
- `SIXEQS_IE_NPHASE`, `Physics.pressure_relaxation = 1`: the reset of the phase internal
  energies after the relaxation left out the term $\alpha_k\rho_kD_k$ of the equation of state,
  so with Mie-Grüneisen and `eref` not zero the phase pressures were not equal afterwards. The
  ideal and stiffened gases give the same results as before.
- `SIXEQS_IE_NPHASE`: the phase sound speed scaled the $B_k$ term of the equation of state by
  $\alpha_k\rho_k$, so every phase with $B_k \neq 0$ (a stiffened gas with $P_{\infty,k} > 0$, or
  Mie-Grüneisen with $p_{\mathrm{ref},k} \neq 0$) had a wrong sound speed, which in SI units made
  the time step too small and the Riemann solvers too diffusive. One formula now serves every
  equation of state, and the mixture uses $\rho c^2 = \sum_k \alpha_k\rho_kc_k^2$ in the
  wave-speed estimates, the time step and the outputs `c` and `Mach`. Runs with ideal gases only
  give the same results as before.
- Viscous and conductive fluxes of `FIVEEQS` and `FIVEEQS_NPHASE` (`-DDIFFUSION=true`): the face
  stencil read cells it does not use, beyond the single ghost cell of
  `FiniteVolume.Scheme = GODUNOV`. These reads are gone. Results are unchanged.
- `run.TimeIntegrator = TVD-RK3` and `RK4`: the stage weights and the flux-register weights of
  the AMR refluxing were truncated decimals. They are now exact quotients. Results with these
  integrators change by round-off; `ForwardEuler` and `TVD-RK2` are unchanged.
- AMR with subcycling, `-DDIFFUSION=true`: a level that a regrid creates during the run got the
  time step $\Delta t_{k-1}/r^2$ but takes only $r$ substeps, so it ended the coarse step behind
  its parent level, which broke conservation and could give NaN in `SIXEQS_IE_NPHASE`. A new level
  now gets $\Delta t_{k-1}/r$, and mass, momentum and energy are conserved when a level appears.
  Builds with `DIFFUSION=false`, runs without subcycling and runs in which no level is created
  are unchanged.
