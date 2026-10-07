# Changelog

Changes that users need to know about, newest first. A change that forces an edit to a case goes
under "Case interface", with a short before-and-after example.

## Unreleased

First public release.

### Added

- Viscous and conductive fluxes for the six-equation models: `-DDIFFUSION=true` now builds
  with `PHYSICS=SIXEQS` and `PHYSICS=SIXEQS_IE_NPHASE`. Each phase has its own stress and its
  own heat flux with its own temperature (Huang, arXiv:2609.18085), the momentum gets the
  mixture stress, and each phase dissipates its own share of the kinetic energy. Mass,
  momentum and total energy stay conservative, also across AMR levels, and an absent phase or a
  massless trace gets nothing. The face viscosity is the harmonic mean of the mixture viscosities
  of the two cells (a negative one counts as 0), so a layered shear flow across a sharp interface
  that lies on a cell face is exact (the arithmetic mean of the five-equation models gives a
  first-order error there; an interface inside a cell keeps that error). With the
  pressure-temperature relaxation of `SIXEQS` the mixture conducts as one medium, with the
  harmonic mean of the cell conductivities, so heat crosses a material interface; with the
  pressure relaxation alone, and always in `SIXEQS_IE_NPHASE` (no temperature relaxation),
  heat does not pass from one phase to another (a known limitation; such a run prints a warning
  when a conductivity is nonzero). The diffusive time step of
  these models uses a true diffusivity. See
  [Viscous and conductive fluxes in the six-equation models](docs/user-guide/models.md#viscous-and-conductive-fluxes-in-the-six-equation-models).
  The previous `SIXEQS_IE_NPHASE` diffusion kernels, which did not compile with the core, are
  replaced. With `DIFFUSION=true`, `SIXEQS` checks that `prob.EOS.cv_1` and `prob.EOS.cv_2` are
  positive.
- Noble-Abel stiffened-gas equation of state (Le Métayer and Saurel 2016), `prob.EOS.EOS = 3`,
  for all four models. The covolume and the heat of formation are new keys, `prob.EOS.b_1`,
  `b_2`, `q_1`, `q_2` for the two-phase models and the lists `prob.EOS.b`, `prob.EOS.q` for the
  N-phase models (default 0), read in `source/Parm.cpp`, so no case file changes. With
  `b = q = 0` the results are identical, bit for bit, to those of the stiffened gas. The covolume
  constraint $1/\rho_k > b_k$ is enforced by a floor on the covolume factor,
  $1 - b_k\rho_k \geq 10^{-6}$. `FiniteVolume.Quad = 1` aborts with EOS 3. New test cases
  `NASG-5Eq`, `NASG-5Eq-N`, `NASG-6Eq`, `NASG-6Eq-N`, and `eos-nasg-b0` variants of the
  `Advection-*` cases. In `SIXEQS_IE_NPHASE` the pressure relaxation keeps the heat of
  formation.
- Test cases (labels of the `Diffusion` directory, see "Changed") `ViscousShockTube-6Eq`,
  `ViscousShockTube-6Eq-N`, `Advection-Viscous-6Eq` and `Advection-Viscous-6Eq-N`, and the
  material-interface tests
  `Couette2Layer-5Eq`, `Couette2Layer-6Eq`, `Couette2Layer-6Eq-N` (two-layer Couette flow, with a
  moving wall through the user boundary condition) and `HeatConduction2Mat-5Eq`,
  `HeatConduction2Mat-6Eq`, `HeatConduction2Mat-6Eq-N` (conduction between two materials), with
  a `mid-cell` variant of `Couette2Layer-6Eq` (the interface inside a cell). `Sod-6Eq` has a
  `water-air-512` variant (the water-air shock tube on 512 uniform cells).
- Test cases `Relaxation-6Eq` and `Relaxation-6Eq-N`: 40 hard states for the pressure relaxation
  (and the pressure-temperature relaxation of `SIXEQS`), one per cell, relaxed once and stepped
  with several equations of state, among them NASG phases next to and past their covolume limit
  and NASG water vapor and liquid water ($q_k \neq 0$ in every phase); the run stops if mass,
  momentum, energy, the volume fraction of a phase that takes no part in the exchange, the sum of
  the volume fractions, the pressure equality or, where they must hold, the energy exchange
  $-P\,d\alpha_k$, the covolume bound $1/\rho_k > b_k$ and, after the pressure-temperature
  relaxation, the temperature equality fail.
- `Physics.pressure_relaxation_trace` (six-equation models, default `0`): a phase with
  $0 < \alpha_k \le$ this value keeps its volume fraction in the pressure relaxation and reaches
  the common pressure through its energy. The closed form of the `SIXEQS` pressure-temperature
  relaxation does not use it, only its fallback to the pressure relaxation. See
  [Safeguarded pressure relaxation](docs/user-guide/models.md#safeguarded-pressure-relaxation).

### Changed

- Test suite (`exec/_Tests`): one directory per problem family. `Advection`, `Diffusion`, `Sod`,
  `NASG`, `Relaxation` and `RayleighTaylor` hold several builds each, one per model
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
  `Advection-5Eq`, which keeps $\alpha_1$ in $[0, 1]$ inside the shapes (the copies of the release
  set $\alpha_1 = -10^{-12}$ there), and write `advection.csv` instead of `interface.csv`, with the
  interface thickness `N_I` of `interface.csv` as its last column (`prob.interface_thickness = 1`
  in their inputs files and in that of `Advection-5Eq`, so the three methods can be compared). `ViscousShockTube-5Eq` and `IsentropicVortex-5Eq` now run WENO5 without
  the quadrature (`FiniteVolume.Quad = 0`), as the six-equation viscous shock tubes do, and a
  variant `quad` of each runs the quadrature WENO5 of their former inputs; the plot folder of
  `IsentropicVortex-5Eq` is `plot/IsentropicVortex`.
  `scripts/test_cases.sh` and `scripts/test_cases.py` build each build of a directory once and
  refuse two builds that would write the same executable; a directory without sections runs as
  before. The suite has 14 directories, 29 builds, 38 base runs and 267 runs.
- AMR without subcycling (`run.do_subcycle = 0`): the levels now take every Runge-Kutta stage
  together, and with `run.do_reflux = 1` the coarse faces covered by a finer level take, at every
  stage and before any level is updated, the area-weighted average of its face fluxes: the
  conservative fluxes with their viscous parts, the face quantities of the
  non-conservative terms (so that a coarse cell's non-conservative terms use the same faces as its
  fluxes), and in the Phase-Field step its fluxes. The coarse-fine ghost cells come from the
  coarser level's stage state (no interpolation in time), and each stage ends with the relaxation
  of the six-equation models and the averaging down. There is no reflux after the step: the update
  is conservative at every stage. Before, each level took all its stages alone and the reflux
  corrected the coarse cells after the step, which with `TVD-RK2`, `TVD-RK3` and `RK4` took the
  volume fraction out of $[0, 1]$ at a sharp interface (up to $9\times10^{-6}$ for a drop crossing
  a coarse-fine boundary at a CFL number of 0.75, the same in every model). Now it stays within
  $[0, 1]$ in every model (`TVD-RK2`, `TVD-RK3` at CFL numbers 0.33 to 0.75 and `RK4`, `GODUNOV`
  and `MUSCL-Minmod`, `pc_interp` and `lincc_interp`, fixed and moving grids, two and three
  levels), with masses and energy conserved to $4\times10^{-15}$, a uniform volume fraction kept
  to $10^{-15}$ under a pressure pulse, and the pressure equilibrium of the relaxed models to
  $3\times10^{-15}$. The error of the isentropic vortex with AMR is 0.2 to 7% smaller with WENO5
  and at most 0.2% larger with `MUSCL-Minmod`. With `run.do_reflux = 0` the stages run the same
  way without the averaged faces. The face fluxes of all levels are held at once during a stage.
  Subcycled runs and runs on one level are unchanged, bit for bit: of the 267 regression runs only
  the two without subcycling and with AMR levels change, `Advection-5Eq@no-subcycle-no-reflux`
  ($2\times10^{-7}$ of $\alpha_1$) and `NonsphericalCollapse-6Eq` ($2\times10^{-4}$ of the largest
  pressure and velocity).

### Fixed

- The configuration log `config_log_<case>.txt` copied `./prob/inputs` whatever inputs file the
  run was given. It now copies the inputs file named on the command line, so a run with
  `prob/inputs.<label>` records its own inputs. Results are unchanged.
- OpenMP (`USE_OMP = TRUE`): in the right-hand side of the hyperbolic and the Phase-Field steps
  each tile multiplied its face fluxes by the face area in place (with `run.do_reflux = 1`), inside
  the parallel loop, while neighboring tiles share their boundary faces (in 3D the default tiles
  are 8 cells long in $y$ and $z$; in 2D a box is one tile unless `fabarray.mfiter_tile_size` makes
  them smaller), so a thread could read a face that another thread had already scaled, or scale it
  twice. Results depended on the thread timing: in a 3D drop test with AMR, 4 threads
  against 1 changed $\alpha_1$ by up to 0.14 after 20 steps, with and without subcycling. The faces
  are now computed once per tile without overlap, the update reads them after all are computed,
  and the area scaling for the flux registers comes after both. The largest wave speed of each
  cell, which sets the time step, is no longer updated by two threads at once (also at
  initialization), and the Phase-Field mechanism fluxes are written once per face. 4 threads now
  give the same bits as 1, and builds without OpenMP are unchanged. Runs with an AMReX linear
  solve (the implicit Phase-Field step, `PhaseField.ID_ExplicitRC = 0`) still depend on the number
  of threads through the order of AMReX's sums, unless `amrex.regtest_reduction = 1`.
- `IsentropicVortex-5Eq` and other cases with `CONVERGENCE` (error norms at the end of the run)
  did not compile with `USE_OMP = TRUE`; the error loop now runs without an OpenMP region. Results
  are unchanged.
- Phase-Field step without subcycling (`PhaseField.phase_field = 1`, `run.do_subcycle = 0`): each
  Phase-Field step also advanced the step counter of the flow, so plotfiles were numbered at twice
  the step number, the last plotfile was written twice (the first one renamed `.old`), and
  `amr.regrid_int` and `amr.chk_int` counted Phase-Field steps. The Phase-Field steps no longer
  count, as with subcycling.
- Six-equation models with a relaxation and AMR (`Physics.pressure_relaxation = 1`, in `SIXEQS`
  also `Physics.pressure_temperature_relaxation = 1`, with `run.do_reflux = 1`): the relaxation
  ran at the end of every Runge-Kutta stage, but the reflux corrected the coarse cells next to a
  finer level after the stages, so those cells ended every step out of equilibrium and started the
  next one from there. With a pressure pulse crossing a coarse-fine boundary,
  $\alpha_1\alpha_2\lvert P_1 - P_2\rvert/P_0$ reached $5\times10^{-3}$ there, against
  $2\times10^{-15}$ elsewhere. The reflux of the hyperbolic step and of the Phase-Field step of
  `SIXEQS` and `SIXEQS_IE_NPHASE` now relaxes again the cells it changed, with the relaxation of the
  stages; the other cells keep their bits, and no volume fraction is clipped. Runs without
  relaxation, without AMR levels or without reflux are unchanged; of the 267 regression runs, the
  91 six-equation runs with a relaxation and AMR change, from round-off amplified in the
  uniform-pressure advection cases to $2\times10^{-4}$ of the largest velocity in
  `ViscousShockTube-6Eq`, where the shock crosses coarse-fine boundaries.
- `SIXEQS`, HLLC Riemann solver: the volume fraction was moved with the velocity of the upwind
  side, $u_K$, and the phase masses with the velocity of the HLLC mass flux,
  $u^*_K = u_K + S_K^\mp\big((S_K - u_K)/(S_K - S^*) - 1\big)$, which is
  $S^*(S_K - u_K)/(S_K - S^*)$ where the star state is used: the star state kept $\alpha_K$
  unscaled (Huang 2026, Eq. 50), so the face velocity of the non-conservative terms, the
  volume-fraction flux at $\alpha = 1$ (Eq. 48), was $u_K$. Where $u_K$ and $S^*$ had opposite
  signs the first-order update of $\alpha_1$ was downwind: $\alpha_1$ left $[0, 1]$ even with
  `GODUNOV` (by up to $10^{-2}$ in a water-air shock tube in a frame moving at 400 m/s, where the
  contact and the fluid next to it start in opposite directions, with NaN when the frame moves at
  600 m/s), a trace gained mass while it lost volume, and in the `NASG-6Eq` tube the covolume floor
  acted ($1 - b\rho$ down to $-0.03$ on 512 cells and $-0.79$ on 1024). The star state now scales
  the volume fraction like the masses, $\alpha^*_K = \alpha_K(S_K - u_K)/(S_K - S^*)$ (Johnsen
  and Colonius 2006, as the other models do), so its flux is $\alpha_K u^*_K$ and the face velocity
  is $u^*_K$, which has the sign of $S^*$: the first-order update of $\alpha_1$ is upwind and keeps
  it within its bounds under the CFL condition, and the `NASG-6Eq` tube no longer reaches the floor
  ($1 - b\rho \geq 0.18$). Interface equilibrium is unchanged ($u^*_K = u$ at a contact in pressure
  and velocity equilibrium). Every `SIXEQS` run with HLLC changes: in the water-air tube of the
  docs by up to 1.7 % of the velocity range in the cells of the shock (L1 errors against the exact
  solution 0.2 to 1.3 % smaller), in the Sod tube by $3\times10^{-6}$ of the density range.
- Pressure relaxation of the six-equation models (`Physics.pressure_relaxation = 1`, and
  `Physics.pressure_temperature_relaxation = 1` in `SIXEQS`): a phase energy that the hyperbolic
  step had left below the minimum of its equation of state ($P_k < -P_{\infty,k}$), a volume
  fraction just outside $[0, 1]$, or a root beyond an asymptote of the relaxation function made it
  fail. In `SIXEQS_IE_NPHASE` the Newton iteration then sat at the pole $P = 0$ of the gas for
  100000 iterations and gave NaN: the water-air shock tube of `Sod-6Eq-N` with any water of
  $\gamma \ge 2$ stopped at the first step, because the air in the cell of the diaphragm had a
  negative energy after the first stage. The iteration could also stall on its absolute test
  $\lvert f\rvert \le 10^{-16}$. Both models keep their own solution wherever it does not fail,
  so results there do not change: finite, in `SIXEQS_IE_NPHASE` a converged root above every
  asymptote of the present phases with the sum of the volume fractions kept (a gas at zero
  pressure lost 70 % of the volume through a $10^{-300}$ guard of the iteration), and in `SIXEQS`
  a root not within $10^{-8}(\lvert P\rvert + \max_k P_{\infty,k})$ of the largest asymptote,
  where the closed form lost its digits (volume fractions off by more than 0.1 for a gas near
  vacuum or a gas trace cavitating in a liquid in tension). A volume fraction outside $[0, 1]$ is
  not a failure: no relaxation clips or rescales a volume fraction, and the bounds come from the
  transport. Elsewhere they use a safeguarded relaxation: phases outside the range of their
  equation of state or without mass keep their volume fraction exactly, and the others share the
  rest of the volume through a Newton iteration with bisection, bracketed between the asymptote
  and a factor 4 of the root of the convex relaxation function, over the whole range of the
  floating-point numbers. The masses and the momentum are unchanged, the total energy is kept,
  $\sum_k\alpha_k = 1$ where a root is found, and no volume fraction is divided by. The
  `SIXEQS_IE_NPHASE` iteration stops after 51 iterations instead of 100000. Where no common
  temperature exists, the pressure-temperature relaxation falls back to the pressure relaxation.
  The closed form of `SIXEQS` and the iteration of `SIXEQS_IE_NPHASE` still give a negative volume
  fraction to a trace whose energy the hyperbolic step left below the minimum of its equation of
  state, and they keep a volume fraction that the hyperbolic step left outside $[0, 1]$.
  `Sod-6Eq-N@water-air` uses the water of the other `Sod` cases again ($\gamma = 4.4$,
  $P_\infty = 6\times10^8$ Pa). `Physics.pressure_relaxation_implicit`, read but never passed to the
  relaxation, is removed. Not fixed: in the water-air shock
  tube, `SIXEQS` with the pressure relaxation alone leaves a cold, over-compressed air layer at the
  contact unless the cell of the diaphragm is mostly water, as on the grid of the docs. On 512
  uniform cells (`Sod-6Eq@water-air-512`, run to $t = 240$ µs) the air at the contact is 18 times
  denser than the shocked air and the shock is 1.1 % of the tube behind. The
  pressure-temperature relaxation and `SIXEQS_IE_NPHASE` do not have it (see
  [Water-air shock tube](docs/verification.md#water-air-shock-tube)).
- `SIXEQS`, `Physics.pressure_temperature_relaxation = 1`: the root of the quadratic of the
  closed form, $(-b + \sqrt{b^2 - 4ac})/(2a)$, lost its digits, since $b > 0$ and $4ac/b^2$ is
  $10^{-3}$ to $10^{-5}$ near 1 bar and $10^{-9}$ or less for a gas near vacuum. The phase
  temperatures after the relaxation differed by $10^{-9}$ to $10^{-8}$ relative near 1 bar and by
  up to 100 % for a gas near vacuum or a gas trace cavitating in a liquid in tension, and with the
  Noble-Abel stiffened gas a liquid could be put at its covolume limit, with a pressure of
  $1.5\times10^9$ Pa instead of about 0. The root is now $-2c/(b + \sqrt{b^2 - 4ac})$ for $b > 0$
  (the temperatures are then equal to round-off), and with the Noble-Abel stiffened gas the closed
  form is kept only where every phase with mass stays below its covolume limit (otherwise the
  pressure relaxation). All runs with this relaxation change by round-off, and more where the
  closed form had lost its digits.
- `SIXEQS_IE_NPHASE`, `Physics.pressure_relaxation = 1`: the reset of the phase internal
  energies after the relaxation left out the term $\alpha_k\rho_kD_k$ of the equation of state,
  so with Mie-Grüneisen and `eref` not zero the phase pressures were not equal afterwards. The
  ideal and stiffened gases ($D_k = 0$) give the same results as before.
- `SIXEQS_IE_NPHASE`: the phase sound speed was computed with $\alpha_k\rho_k$ in place of
  $\alpha_k$ and $\alpha_kB_k$ in place of $B_k$, which gave
  $\alpha_k\rho_kc_k^2 = ((A_k+1)\alpha_kP_k + \alpha_k\rho_k\,\alpha_kB_k)/A_k$. Every phase with
  $B_k \neq 0$ (a stiffened gas with $P_{\infty,k} > 0$, or Mie-Grüneisen with
  $p_{\mathrm{ref},k} \neq 0$) had its $B_k$ term scaled by $\alpha_k\rho_k$. In SI units the sound
  speed of nearly pure water came out about $\sqrt{\rho_k} \approx 30$ times too large, so the time
  step was that much smaller and the Riemann solvers more diffusive. One formula now serves every equation of state,
  $\alpha_k\rho_kc_k^2 = ((A_k+1)\alpha_kP_k + \alpha_kB_k)/(A_k\xi_k)$ with $\xi_k = 1$ except for
  the Noble-Abel stiffened gas, and the mixture uses $\rho c^2 = \sum_k \alpha_k\rho_kc_k^2$ in the
  wave-speed estimates of the HLLC, HLL and LLF solvers, the time step and the outputs `c` and
  `Mach`. Runs with $B_k = 0$ for every phase (ideal gases) give the same results as before.
- Viscous and conductive fluxes of `FIVEEQS` and `FIVEEQS_NPHASE` (`-DDIFFUSION=true`): the face
  stencil also read the cells one cell beyond the face cell along the face normal, which the
  face does not use; on the last face of a box they lie two cells outside it, beyond the single
  ghost cell of `FiniteVolume.Scheme = GODUNOV`. These reads are gone, as in the six-equation
  models. Results are unchanged.
- `run.TimeIntegrator = TVD-RK3` and `RK4`: the stage weights and the flux-register weights of
  the AMR refluxing were truncated decimals (`0.166666666667`, `0.333333333333`,
  `0.666666666667`), off by up to $2\times10^{-12}$ relative. They are now exact quotients such as
  `Real(2.0)/Real(3.0)`. Results with these integrators change by round-off-sized perturbations;
  `ForwardEuler` and `TVD-RK2` are unchanged.
- AMR with subcycling, `-DDIFFUSION=true`: a level that a regrid creates during the run got the
  time step $\Delta t_{k-1}/r^2$ but takes only $r$ substeps, so it ended the coarse step behind
  its parent level. Its next fills then extrapolated in time, which changed the conserved totals
  (a mass change of order $10^{-4}$ in the coarse step after the level appeared) and could give
  NaN in `SIXEQS_IE_NPHASE`. A new level now gets $\Delta t_{k-1}/r$, as everywhere else in the
  time-step code, and mass, momentum and energy are conserved to round-off when a level appears.
  Builds with `DIFFUSION=false`, runs without subcycling and runs in which no level is created
  are unchanged.
