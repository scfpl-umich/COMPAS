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
  momentum and total energy stay conservative, also across AMR levels, and an absent phase
  gets nothing. The diffusive time step of these models uses a true diffusivity. See
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
- Test cases `ShearDecay-6Eq`, `ShearDecay-6Eq-N`, `HeatConduction-6Eq`,
  `HeatConduction-6Eq-N`, `ViscousShockTube-6Eq`, `ViscousShockTube-6Eq-N`,
  `Advection-Viscous-6Eq` and `Advection-Viscous-6Eq-N`.

### Fixed

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
