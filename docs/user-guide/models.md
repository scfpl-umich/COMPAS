# Models and equations

COMPAS solves diffuse-interface models of compressible multiphase flow. The model and the viscous
and conductive fluxes are chosen when the case is compiled, with `-DPHYSICS=` and `-DDIFFUSION=`
in the `GNUmakefile` ([Compile-time options](compile-options.md)). The other terms below and the
numerical methods are switched at run time from the inputs file ([Run-time options](inputs.md)).

## Notation and equations of state

Phase $k = 1,\dots,N$ has a volume fraction $\alpha_k$, density $\rho_k$, specific
internal energy $e_k$ and pressure $P_k$, with $\sum_k \alpha_k = 1$. All phases move
with one velocity $\vec u$, in every model. The mixture quantities are

$$
\rho = \sum_k \alpha_k\rho_k, \qquad
Y_k = \frac{\alpha_k\rho_k}{\rho}, \qquad
\rho e = \sum_k \alpha_k\rho_k e_k, \qquad
\rho E = \rho e + \tfrac12\rho\lvert\vec u\rvert^2 .
$$

Interfaces are diffuse: $\alpha_k$ changes smoothly across a few cells, and no interface
tracking is needed.

Every phase has an equation of state of the form

$$
\rho_k e_k = \left(1 - C_k\rho_k\right)\left(A_k P_k + B_k\right) + \rho_k D_k ,
$$

which covers the four equations of state in COMPAS. The run-time parameter `prob.EOS.EOS`
selects one for all phases.

| `prob.EOS.EOS` | $A_k$ | $B_k$ | $C_k$ | $D_k$ |
|---|---|---|---|---|
| `0`, ideal gas | $1/(\gamma_k-1)$ | $0$ | $0$ | $0$ |
| `1`, stiffened gas | $1/(\gamma_k-1)$ | $\gamma_k P_{\infty,k}/(\gamma_k-1)$ | $0$ | $0$ |
| `2`, Mie–Grüneisen | $1/\Gamma_k$ | $-p_{\mathrm{ref},k}/\Gamma_k$ | $0$ | $e_{\mathrm{ref},k}$ |
| `3`, Noble-Abel stiffened gas | $1/(\gamma_k-1)$ | $\gamma_k P_{\infty,k}/(\gamma_k-1)$ | $b_k$ | $q_k$ |

The constants $\gamma_k$, $P_{\infty,k}$, $\Gamma_k$, $p_{\mathrm{ref},k}$,
$e_{\mathrm{ref},k}$, $b_k$, $q_k$ and $c_{v,k}$ are the inputs `gamma`, `pinf`, `gGamma`, `pref`,
`eref`, `b`, `q` and `cv`, listed under
[Equation of state and materials](inputs.md#equation-of-state-and-materials).

Where a temperature is needed, it comes from

$$
\rho_k e_k = \rho_k c_{v,k} T_k + \left(1 - C_k\rho_k\right)\frac{B_k}{1+A_k} + \rho_k D_k .
$$

### Noble-Abel stiffened gas

The Noble-Abel stiffened gas (NASG, Le Métayer and Saurel 2016) adds to the stiffened gas the
covolume $b_k$, the volume of the molecules, and the heat of formation $q_k$:

$$
\rho_k e_k = \frac{\left(1 - b_k\rho_k\right)\left(P_k + \gamma_k P_{\infty,k}\right)}{\gamma_k - 1} + \rho_k q_k ,
\qquad
\rho_k c_k^2 = \frac{\gamma_k\left(P_k + P_{\infty,k}\right)}{1 - b_k\rho_k} ,
\qquad
T_k = \frac{\left(P_k + P_{\infty,k}\right)\left(1/\rho_k - b_k\right)}{(\gamma_k - 1)\,c_{v,k}} .
$$

With $b_k = q_k = 0$ it is the stiffened gas, and the results are then identical, bit for bit,
to those of `prob.EOS.EOS = 1`. In the mixture equations below the covolume enters through
$\phi_k = \alpha_k - b_k\alpha_k\rho_k$, which replaces $\alpha_k$ wherever it multiplies
$A_k$ or $B_k$, so that $\alpha_k\rho_ke_k = \phi_k(A_kP_k + B_k) + \alpha_k\rho_k q_k$.

The equation of state needs $1/\rho_k > b_k$, that is $\xi_k = 1 - b_k\rho_k > 0$. Numerical
states can break it, mostly in the traces of a phase at interfaces, where the mass and the volume
fraction of the phase are not transported in proportion. COMPAS bounds the covolume factor
where it is formed, with a floor: $\xi_k \geq 10^{-6}$ in the six-equation models and
$\phi_k \geq 10^{-6}\,\alpha_k$ in the five-equation models (`NASG_Xi_Min` in
`include/EquationOfState_NASG.H`). The floor is applied only to phases with $b_k > 0$. It keeps
the state usable and leaves the conservative variables untouched, but where it acts the equation
of state is no longer linear in $(\alpha_k, \alpha_k\rho_k)$ at fixed pressure, so pressure
equilibrium at interfaces is no longer kept to round-off in those cells. In the six-equation
models the pressure relaxation computes the relaxed volume fractions from the unbounded
$\phi_k$, which moves such a phase back toward $1/\rho_k > b_k$. The HLLC solvers move the volume
fraction and the mass of a phase with the same velocity, so a trace does not gain mass while it
loses volume: in the air-water shock tube of `NASG-6Eq`, next to its strong pressure jump at the
interface, the smallest unbounded $\xi_k$ is 0.18 on 512 and on 1024 cells, and the floor does
not act. Where the floor acts, conservation is unaffected, but the sound speed of a
floored phase, $c_k^2 \propto 1/\xi_k$, shortens the time step. The bound-preserving limiter (`FiniteVolume.ID_Bound = 1`) works with NASG,
but it bounds the energy, not $\phi_k$.

**Reference.** O. Le Métayer and R. Saurel, The Noble-Abel stiffened-gas equation of state,
*Phys. Fluids* 28 (2016) 046102, [doi:10.1063/1.4945981](https://doi.org/10.1063/1.4945981).

## Five-equation model: `FIVEEQS` and `FIVEEQS_NPHASE`

The phases share one pressure, $P_k = P$.

$$
\begin{aligned}
\frac{\partial (\alpha_k\rho_k)}{\partial t} + \nabla\cdot\left(\alpha_k\rho_k \vec{u}\right) &= 0\\
\frac{\partial (\rho \vec{u})}{\partial t} + \nabla\cdot\left(\rho \vec{u}\otimes\vec{u} + P\,\mathbb{I}\right) &= \nabla\cdot\boldsymbol{\tau}\\
\frac{\partial (\rho E)}{\partial t} + \nabla\cdot\left[(\rho E+P)\,\vec{u}\right] &= \nabla\cdot\left(\boldsymbol{\tau}\cdot\vec{u} - \vec{q}\right)\\
\frac{\partial \alpha_k}{\partial t} + \vec{u}\cdot\nabla\alpha_k &= \alpha_k K_k\,\nabla\cdot\vec{u}
\end{aligned}
$$

with $k = 1,\dots,N$ for the masses and $k = 1,\dots,N-1$ for the volume fractions. The
mixture equation of state gives the pressure,

$$
P = \frac{\rho e - \sum_k \left(\phi_k B_k + \alpha_k\rho_k D_k\right)}{\sum_k \phi_k A_k} ,
\qquad \phi_k = \alpha_k\left(1 - C_k\rho_k\right) ,
$$

which needs no division by $\alpha_k$; $\phi_k = \alpha_k$ except with the Noble-Abel stiffened
gas.

`FIVEEQS` has two phases. `FIVEEQS_NPHASE` has $N$ phases, set with `-DNPHASE=`$N$ in the
`GNUmakefile`.

`Physics.source_term` sets the volume-fraction source.

- `0` (default): $K_k = 0$, so $\alpha_k$ is only advected (Allaire et al. 2002).
- `1`: $K_k = \dfrac{\rho c^2}{\rho_k c_k^2} - 1$, with
  $\dfrac{1}{\rho c^2} = \sum_j \dfrac{\alpha_j}{\rho_j c_j^2}$ (Kapila et al. 2001).
  With the Noble-Abel stiffened gas this is computed as
  $\alpha_k(1 + K_k) = \phi_k\,\rho c^2/(\gamma_k(P + P_{\infty,k}))$ and
  $1/(\rho c^2) = \sum_j \phi_j/(\gamma_j(P + P_{\infty,j}))$, again without a division by
  $\alpha_k$.

With `1`, each phase in a mixture cell is compressed in inverse proportion to its stiffness
$\rho_k c_k^2$, so all phases feel the same pressure change; use it when shocks cross interfaces
between different materials.

With `-DDIFFUSION=true` the viscous stress and heat flux are

$$
\boldsymbol\tau = \mu\left(\nabla\vec u + \nabla\vec u^{T}\right) + \left(\mu_B - \tfrac23\mu\right)(\nabla\cdot\vec{u})\,\mathbb{I},
\qquad
\vec q = -\kappa\nabla T,
$$

with $\mu = \sum_k\alpha_k\mu_k$, and $\mu_B$ and $\kappa$ mixed the same way. $T$ comes
from the mixture equation of state in its temperature form. Without it,
$\boldsymbol\tau = 0$ and $\vec q = 0$. The phase values $\mu_k$, $\mu_{B,k}$ and
$\kappa_k$ are the inputs `mu`, `muB` and `kappa`.

**References.** Two-phase model with Phase-Field: Z. Huang and E. Johnsen, *J. Comput. Phys.*
488 (2023) 112195, [doi:10.1016/j.jcp.2023.112195](https://doi.org/10.1016/j.jcp.2023.112195).
N-phase model with Phase-Field: Z. Huang and E. Johnsen, *J. Comput. Phys.* 501 (2024) 112801,
[doi:10.1016/j.jcp.2024.112801](https://doi.org/10.1016/j.jcp.2024.112801).

## Six-equation model: `SIXEQS`

Two phases, each with its own pressure $P_k$ and total energy
$E_k = e_k + \tfrac12\lvert\vec u\rvert^2$ (Pelanti and Shyue 2014).

$$
\begin{aligned}
\frac{\partial \alpha_1}{\partial t} + \vec{u}\cdot\nabla\alpha_1 &= \mathcal{R}\\
\frac{\partial (\alpha_k\rho_k)}{\partial t} + \nabla\cdot\left(\alpha_k\rho_k \vec{u}\right) &= 0\\
\frac{\partial (\rho \vec{u})}{\partial t} + \nabla\cdot\left(\rho \vec{u}\otimes\vec{u} + P\,\mathbb{I}\right) &= 0\\
\frac{\partial (\alpha_1\rho_1 E_1)}{\partial t} + \nabla\cdot\left[\alpha_1(\rho_1E_1 + P_1)\,\vec{u}\right] &= \phantom{-}\Sigma - P_I\,\mathcal{R}\\
\frac{\partial (\alpha_2\rho_2 E_2)}{\partial t} + \nabla\cdot\left[\alpha_2(\rho_2E_2 + P_2)\,\vec{u}\right] &= -\Sigma + P_I\,\mathcal{R}
\end{aligned}
$$

with $P = \alpha_1P_1 + \alpha_2P_2$ and the non-conservative exchange term

$$
\Sigma = \vec{u}\cdot\left[\,Y_2\nabla(\alpha_1P_1) - Y_1\nabla(\alpha_2P_2)\right] .
$$

- Each phase has its own equation of state,
  $\alpha_k\rho_ke_k = \xi_k A_k\,\alpha_kP_k + \xi_k\alpha_kB_k + \alpha_k\rho_kD_k$, with
  $\xi_k = 1 - C_k\rho_k$ ($\xi_k = 1$ except with the Noble-Abel stiffened gas).
- The two energy equations add up to the conservative equation for $\rho E$, so total energy
  is conserved.
- The hyperbolic step follows Huang (2026): the face velocity $\hat u$ of the non-conservative
  terms ($\alpha_1\nabla\cdot\hat u$ and $\Sigma$) is the numerical flux of $\alpha_1$ evaluated
  at $\alpha_1 = 1$, and the products $\alpha_kP_k\vec u$ in $\Sigma$ are taken from the energy
  fluxes. In the HLLC solver (`Physics.RiemannSolver = "HLLC"`) the star states scale the volume
  fraction like the partial densities, $\alpha_K^* = \alpha_K(S_K - u_K)/(S_K - S^*)$ for
  $K = L, R$ (Johnsen and Colonius 2006, as in the five-equation models; Huang keeps
  $\alpha_K^* = \alpha_K$). The volume fraction then moves with the velocity of the mass flux,
  $\hat u = u_K^* = S^*(S_K - u_K)/(S_K - S^*)$ where the star state is used and $u_K$ where the
  flux is $F(U_K)$, which has the sign of the contact speed $S^*$, so the first-order update of
  $\alpha_1$ is upwind and keeps it within its bounds under the CFL condition. A contact in pressure
  and velocity equilibrium gives $S^* = u_K^* = u$, so interface equilibrium is kept. In the HLL and
  LLF solvers $\hat u$ is likewise the volume-fraction flux at $\alpha_1 = 1$.

![Bounds of the volume fraction with the two HLLC star states](../media/figures/hllc_alpha_bounds.png)

*Volume fraction of the `SIXEQS` HLLC solver in the
[water-air shock tube](../verification.md#water-air-shock-tube) seen from a frame moving at
400 m/s (both sides at $u_0 = -400$ m/s, so the contact moves against the fluid next to it;
256 cells, no relaxation, each material a trace of $10^{-6}$ on the other side), with the star
state of Huang, $\alpha_K^* = \alpha_K$, and that of COMPAS, which scales $\alpha$ like the
partial densities. (a) GODUNOV and (b) MUSCL: the smallest $\alpha_1$ and the largest
$\alpha_1 - 1$ against the step. (c) $\alpha_1$ and $\alpha_1\rho_1$ next to the contact after
three GODUNOV steps. With the scaled star state $\alpha_1$ stays within its initial bounds,
$[10^{-6}, 1 - 10^{-6}]$.*

The HLLC flux moves the partial densities with $u_K^*$, so $\alpha$ must move with the same face
velocity: moved with $u_K$, it is updated from the downwind cell wherever $u_K$ and $S^*$ differ in
sign, here falling to $-5.9 \times 10^{-3}$ while $\alpha_1\rho_1$ stays positive.

- The relaxation $\mathcal{R} = \mu_P(P_1 - P_2)$ is instantaneous, $\mu_P \to \infty$, and
  is applied after every Runge–Kutta stage, and with AMR and subcycling again in the coarse cells
  that the reflux corrects after the stages, so that a step ends in equilibrium there too. Without
  subcycling the coarse cells next to a finer level get its fluxes within each stage
  ([AMR time stepping](inputs.md#time-stepping)), so they are relaxed with all the others.
- `Physics.pressure_relaxation = 1` relaxes to $P_1 = P_2$, with $P_I$ the relaxed pressure.
- `Physics.pressure_temperature_relaxation = 1` relaxes to $P_1 = P_2$ and $T_1 = T_2$.
- Both relaxations have closed forms, also with the Noble-Abel stiffened gas: a quadratic
  for the relaxed pressure, written with $\phi_k = \xi_k\alpha_k$ and the phase energies
  without the heat of formation, which reduces to the stiffened-gas one when $b_k = 0$. The
  pressure-temperature relaxation takes the larger root of $ap^2 + bp + c = 0$ as
  $-2c/(b + \sqrt{b^2 - 4ac})$ when $b > 0$, the usual case: $(-b + \sqrt{b^2 - 4ac})/(2a)$ loses
  the digits of a pressure that is small against $P_{\infty,k}$.
- The closed form is used wherever it does not fail: a finite result (and not the spurious root
  $P = -P_{\infty,k}$ of an absent phase), with its root farther than
  $10^{-8}(\lvert P\rvert + \max_k P_{\infty,k})$ from the largest asymptote $-P_{\infty,k}$ of the
  present phases. Closer, the quadratic loses most of its digits (a gas near vacuum, or a gas trace
  cavitating in a liquid in tension: volume fractions off by more than 0.1). There the pressure
  relaxation is the [safeguarded relaxation](#safeguarded-pressure-relaxation). The
  pressure-temperature relaxation, where no common temperature exists (a mixture energy too low
  for $T > 0$, or with the Noble-Abel stiffened gas no volume left for a phase above its covolume
  $b_k\alpha_k\rho_k$), falls back to the pressure relaxation.
- With neither, the phases keep their own pressures. Both are off by default.
- No relaxation clips or rescales a volume fraction: a volume fraction outside $[0, 1]$ is not a
  failure of the closed form, which keeps it (and a trace whose energy the hyperbolic step left
  below the minimum of its equation of state can get a negative volume fraction from it). The
  bounds come from the transport: the first-order update keeps $\alpha_1$ within the values of its
  neighbors (see the HLLC solver above), and so does MUSCL with a TVD limiter, while WENO5 (in
  `Advection-6Eq` down to $\alpha_1 = -2\times10^{-3}$ in 8 steps) and to a far smaller extent
  WENO5B can take it slightly outside $[0, 1]$, as can the refluxing of a subcycled step
  (`run.do_subcycle = 1`, `run.do_reflux = 1`) with the multistage integrators `TVD-RK2`,
  `TVD-RK3` and `RK4`: $2.9\times10^{-6}$ in the first step of `Advection-6Eq` with `GODUNOV` and
  `TVD-RK3`, none without the reflux, with either interpolater, and none seen with `ForwardEuler`.
  It grows with the time step: for a sharp drop advected across a coarse-fine boundary with
  `TVD-RK3`, none at a CFL number $(|u| + c)\Delta t/\Delta x$ of 0.17, $2\times10^{-6}$ at 0.33 and
  $10^{-4}$ at 0.67, in the first steps only. This is the same in the five-equation models, to the
  last digit: the reflux replaces the coarse flux at the coarse-fine face after the step, but not
  the stage values of the coarse cell, which were built with the coarse inflow and from which its
  outflow (through its other faces, and into a finer level through the ghost cells filled from it)
  was computed. Where the fine inflow is much smaller than the coarse one (the steep tail of a sharp
  interface), the corrected value is not a convex combination of the old values. The reflux keeps
  the masses and the energy conserved, the volume fractions summing to 1, and a uniform pressure
  and velocity, to round-off. Without subcycling (`run.do_subcycle = 0`) there is no reflux: the
  coarse cells take the fluxes of the finer level within every stage
  ([AMR time stepping](inputs.md#time-stepping)), and the same drop stays within $[0, 1]$ in
  every model, with `TVD-RK2`, `TVD-RK3` up to a CFL number of 0.75 and `RK4`, where the reflux
  after the step gave up to $9\times10^{-6}$. Conservation and the sum of the volume fractions are
  kept to round-off, and a uniform pressure and velocity as well as on one level at the same CFL
  number.
- With `-DDIFFUSION=true` the model has viscous and conductive fluxes, described under
  [Viscous and conductive fluxes in the six-equation models](#viscous-and-conductive-fluxes-in-the-six-equation-models).

**References.** Z. Huang, A consistent and conservative Phase-Field method for compressible
multiphase flows with the six-equation model, arXiv:2609.18085 (2026),
[doi:10.48550/arXiv.2609.18085](https://doi.org/10.48550/arXiv.2609.18085). HLLC star state of the
volume fraction: E. Johnsen and T. Colonius, Implementation of WENO schemes in compressible
multicomponent flow problems, *J. Comput. Phys.* 219 (2006) 715-732,
[doi:10.1016/j.jcp.2006.04.018](https://doi.org/10.1016/j.jcp.2006.04.018).

## N-phase six-equation model: `SIXEQS_IE_NPHASE`

$N$ phases with their own pressures, written with the phase internal energies and the
mixture total energy (Saurel et al. 2009).

$$
\begin{aligned}
\frac{\partial \alpha_k}{\partial t} + \vec{u}\cdot\nabla\alpha_k &= \mathcal{R}_k\\
\frac{\partial (\alpha_k\rho_k)}{\partial t} + \nabla\cdot\left(\alpha_k\rho_k \vec{u}\right) &= 0\\
\frac{\partial (\alpha_k\rho_k e_k)}{\partial t} + \nabla\cdot\left(\alpha_k\rho_k e_k\, \vec{u}\right) + \alpha_kP_k\,\nabla\cdot\vec{u} &= -P_I\,\mathcal{R}_k\\
\frac{\partial (\rho \vec{u})}{\partial t} + \nabla\cdot\left(\rho \vec{u}\otimes\vec{u} + P\,\mathbb{I}\right) &= 0\\
\frac{\partial (\rho E)}{\partial t} + \nabla\cdot\left[(\rho E+P)\,\vec{u}\right] &= 0
\end{aligned}
$$

for $k = 1,\dots,N$, with $P = \sum_k \alpha_kP_k$ and $\sum_k\mathcal{R}_k = 0$.

- All $N$ volume fractions, partial densities and internal energies are evolved, plus one
  momentum and one total energy.
- The internal energy equations are non-conservative. The total energy equation is
  conservative, and the relaxation makes the two agree.
- The wave-speed estimates of the Riemann solvers, the time step and the output `c` use the
  frozen sound speed $\rho c^2 = \sum_k \alpha_k\rho_kc_k^2$, with
  $\rho_kc_k^2 = ((A_k+1)P_k + B_k)/(A_k\xi_k)$ for every equation of state, which is
  $\gamma_k(P_k + P_{\infty,k})/\xi_k$ for the stiffened gas and the Noble-Abel stiffened gas.
- `Physics.pressure_relaxation = 1`: after every stage (and with AMR and subcycling again in the
  coarse cells that the reflux corrects) a Newton iteration finds the common pressure $P$ from

  $$
  \sum_k \frac{\alpha_k\,(P_k - P)}{\rho_k c_k^2(P)} = 0 ,
  $$

  the volume fractions follow with $P_I = P$, and the phase internal energies are reset from
  $\rho e = \rho E - \tfrac12\rho\lvert\vec u\rvert^2$. It is off by default. With the
  Noble-Abel stiffened gas, $\rho_kc_k^2(P) = \gamma_k(P + P_{\infty,k})/\xi_k$, so the equation
  reads $\sum_k \phi_k(P_k - P)/(\gamma_k(P + P_{\infty,k})) = 0$, and the reset keeps the
  heat-of-formation terms $\alpha_k\rho_kq_k$.
  The Newton iteration starts from $\lvert\sum_k\alpha_kP_k\rvert$, halved until the function is
  positive, and stops at $\lvert f\rvert \le 10^{-16}$, where $\lvert f'\rvert < 10^{-14}$, or after 51
  iterations. Its result is used
  where it does not fail: $\lvert f\rvert \le 10^{-12}$ at a pressure above every $-P_{\infty,k}$ of
  the present phases, and relaxed volume fractions whose sum is the sum before to
  $10^{-12}$ ($f$ has a $10^{-300}$ guard in its denominators that the volume fractions do not, so
  a gas at zero pressure could pass the test on $f$ and lose 70 % of the volume). Elsewhere (for
  instance a phase energy below the minimum of its equation of state, which can leave the function
  without a root above the asymptotes) the volume fractions come from the
  [safeguarded relaxation](#safeguarded-pressure-relaxation); the reset is the same. As in
  `SIXEQS`, a volume fraction outside $[0, 1]$ is kept, not clipped.
- $N$ is set with `-DNPHASE=`$N$ in the `GNUmakefile`.
- With `-DDIFFUSION=true` the model has viscous and conductive fluxes, described in the next
  section.

**Reference.** The same paper as `SIXEQS`: Z. Huang, arXiv:2609.18085 (2026),
[doi:10.48550/arXiv.2609.18085](https://doi.org/10.48550/arXiv.2609.18085).

## Safeguarded pressure relaxation

Both six-equation models fall back on it where their own solution of the pressure relaxation fails
(see above). With the phase energies exchanged by $-P\,d\alpha_k$ at the final pressure
$P$, the relaxed volume fractions are (Huang 2026, Eq. 32)

$$
\alpha_k(P) = \alpha_k + \frac{\phi_k\,(P_k - P)}{\gamma_k\,(P + P_{\infty,k})},
$$

with $\phi_k = \xi_k\alpha_k$, and $\gamma_k = 1 + \Gamma_k$,
$P_{\infty,k} = -p_{\mathrm{ref},k}/(1 + \Gamma_k)$ for the Mie-Grüneisen equation of state. The
products $\phi_kP_k$ come from the phase energies, so no volume fraction is divided by.

- A phase takes part in the exchange when $\alpha_k > 0$, $\alpha_k\rho_k > 0$,
  $\phi_k(P_k + P_{\infty,k}) > 0$ (that is $P_k > -P_{\infty,k}$, except for a NASG phase past
  its covolume limit, $\phi_k < 0$) and $\alpha_k$ is above `Physics.pressure_relaxation_trace`
  (default 0; the closed form of the `SIXEQS` pressure-temperature relaxation does not use the
  threshold, only its fallback to the pressure relaxation does).
  The others (absent, a volume fraction that the hyperbolic step left negative, a phase without
  mass, an energy below the minimum of the equation of state, or a trace below the threshold) keep
  their volume fraction exactly and reach the common pressure through their energy.
- The phases that take part share the rest of the volume, $1 - \sum$ of the others. Their
  relaxation function is convex and decreasing above the largest $-P_{\infty,k}$ among them, where
  it has exactly one root (Huang 2026, Theorems 3.1 and 3.2), at which each of them has a positive
  volume fraction. The root is bracketed by factors of 4 from the mean pressure of
  these phases, over the whole range of the floating-point numbers (a root next to the asymptote
  of a $10^{-200}$ trace included), and found by Newton iterations, with a bisection where an
  iterate leaves the bracket, to a relative tolerance of $10^{-14}$ in $P + P_{\infty,k}$, or until
  the function is zero to its round-off (at most 100 iterations).
- The phase energies are then reset at the pressure given by the total energy: $\rho E$ in
  `SIXEQS_IE_NPHASE`, $E_1 + E_2$ in `SIXEQS`. The masses and the momentum are unchanged, the
  total energy is kept to round-off, and $\sum_k\alpha_k = 1$. No volume fraction is clipped or
  rescaled: where no phase takes part, or the function has no root (an input whose volume fractions
  add up to far more than 1), every volume fraction is kept.
- If the total energy is below the energy that the phases need at the relaxed volume fractions,
  the common pressure falls below $-P_{\infty,k}$ of a phase. Energy conservation leaves no other
  choice; the hyperbolic step should not produce such a state.

The two base runs of `exec/_Tests/Relaxation`, `Relaxation-6Eq` and `Relaxation-6Eq-N`, relax 40
hard states with each model (traces compressed and in tension, absent phases, negative volume
fractions, energies at and below the minimum of the equation of state, pressure ratios up to $10^7$,
NASG phases next to and past the covolume limit, the two cells of the water-air shock tube where the
N-phase iteration used to fail, a gas at zero energy, and gas traces cavitating next to the gas
asymptote) and check conservation, that a phase taking no part keeps its volume fraction to the bit,
that the volume fractions add up to 1, and pressure equality, and, on the states where they must
hold, the energy exchange $-P\,d\alpha_k$, a relaxed pressure above every $-P_{\infty,k}$, the
covolume bound $1/\rho_k > b_k$ and, after the pressure-temperature relaxation, equal phase
temperatures. The pressure equality alone holds for any volume fractions, since the phase energies
are reset at one pressure.

## Viscous and conductive fluxes in the six-equation models

With `-DDIFFUSION=true`, `SIXEQS` and `SIXEQS_IE_NPHASE` give each phase its own viscous stress
and heat flux (Huang 2026),

$$
\boldsymbol\tau_k = \mu_k\left(\nabla\vec u + \nabla\vec u^{T}\right) + \left(\mu_{B,k} - \tfrac23\mu_k\right)(\nabla\cdot\vec{u})\,\mathbb{I},
\qquad
\vec q_k = -\kappa_k\nabla T_k ,
$$

with $\mu_{B,k}$ the bulk viscosity, as in the five-equation models (the coefficient of
$\nabla\cdot\vec u$ in the paper is $\mu_{B,k} - \tfrac23\mu_k$), and with $T_k$ the temperature of
phase $k$ from its own equation of state. The momentum gets the mixture stress
$\boldsymbol\tau = \sum_k\alpha_k\boldsymbol\tau_k$, and the energies get, in `SIXEQS` (the
two-phase model) and in `SIXEQS_IE_NPHASE` (the N-phase model),

$$
\begin{aligned}
\frac{\partial (\rho \vec{u})}{\partial t} + \dots &= \nabla\cdot\boldsymbol\tau\\
\frac{\partial (\alpha_k\rho_k E_k)}{\partial t} + \dots &= Y_k\,\vec u\cdot(\nabla\cdot\boldsymbol\tau) + \alpha_k\boldsymbol\tau_k : \nabla\vec u - \nabla\cdot(\alpha_k\vec q_k) && \text{(two-phase model)}\\
\frac{\partial (\alpha_k\rho_k e_k)}{\partial t} + \dots &= \alpha_k\boldsymbol\tau_k : \nabla\vec u - \nabla\cdot(\alpha_k\vec q_k) && \text{(N-phase model)}\\
\frac{\partial (\rho E)}{\partial t} + \dots &= \nabla\cdot\Big(\boldsymbol\tau\cdot\vec u - \sum_k\alpha_k\vec q_k\Big) && \text{(N-phase model)}
\end{aligned}
$$

Each phase dissipates $\alpha_k\boldsymbol\tau_k : \nabla\vec u$ with its own viscosity, and in
`SIXEQS` the work of the mixture viscous force is shared by mass fraction. The discretization:

- **Momentum.** The face flux $\boldsymbol\tau_f\cdot\vec e_d$ uses the face velocity gradient of
  the five-equation models and the face viscosities, the harmonic means
  $2\mu_m\mu_c/(\mu_m + \mu_c)$ of the mixture viscosities $\mu = \sum_k\alpha_k\mu_k$ of the two
  cells (and likewise for $\mu_B = \sum_k\alpha_k\mu_{B,k}$); see
  [Material interfaces](#material-interfaces).
- **Viscous energy.** The conservative face flux $\vec u_f\cdot\boldsymbol\tau_f\cdot\vec e_d$ goes
  to $\rho E$ in `SIXEQS_IE_NPHASE`, and each phase internal energy gets the cell source
  $\alpha_k\Phi_k$, with
  $\Phi_k = \mu_k(\partial_iu_j + \partial_ju_i)\partial_ju_i + (\mu_{B,k} - \tfrac23\mu_k)(\nabla\cdot\vec u)^2 \ge 0$
  from central differences. `SIXEQS` stores no mixture energy: the whole flux goes to the energy
  of phase 1, and phase 2 then takes
  $Y_2[\nabla\cdot(\vec u\cdot\boldsymbol\tau)]_h + \alpha_2\Phi_2 - Y_2\Phi$ from it, with
  $\Phi = \sum_k\alpha_k\Phi_k$ and $[\nabla\cdot(\vec u\cdot\boldsymbol\tau)]_h$ the divergence of
  the same face fluxes. Phase $k$ thus gets
  $Y_k[\nabla\cdot(\vec u\cdot\boldsymbol\tau)]_h + \alpha_k\Phi_k - Y_k\Phi$. In these shares
  $\alpha_k$ counts only where the phase has mass ($\alpha_k > 0$ and $\alpha_k\rho_k > 0$), so a
  massless trace left by round-off gets no energy (it would otherwise get $\alpha_k\Phi_k$, of
  either sign with the sign of the round-off).
- **Heat, phase by phase.** Without the pressure-temperature relaxation each phase conducts
  through itself, with the face flux $w_k\,\kappa_k\,(T_{k,c} - T_{k,m})/\Delta x$ into its
  energy, $w_k$ the harmonic mean $2\alpha_{k,m}\alpha_{k,c}/(\alpha_{k,m} + \alpha_{k,c})$ of the
  volume fractions of the two cells. $T_k$ is used only where the phase is present,
  $\alpha_k \ge 10^{-10}$ and $\alpha_k\rho_k > 0$, and $w_k = 0$ unless it is present in both
  cells. In `SIXEQS_IE_NPHASE` the sum of the phase heat fluxes also goes to $\rho E$.
- **Heat with the pressure-temperature relaxation** (`SIXEQS`,
  `Physics.pressure_temperature_relaxation = 1`). The relaxation gives the phases of a cell one
  temperature, and the relaxed state depends only on the phase masses, the momentum and the
  total energy, not on how the energy is shared between the phases. The mixture then conducts as
  one medium, $\vec q = -\kappa\nabla T$ with $\kappa = \sum_k\alpha_k\kappa_k$ (the sum of the
  phase heat fluxes): the face flux is $\kappa_f(T_c - T_m)/\Delta x$, with $\kappa_f$ the
  harmonic mean of the conductivities of the two cells and $T$ the mean of the temperatures of
  the phases present weighted by their heat capacities $\alpha_k\rho_kc_{v,k}$, the temperature
  they reach by exchanging heat at fixed volumes (after the relaxation the phase temperatures are
  equal, so the weights matter only in cells where the relaxation fell back to the pressure
  relaxation and in coarse-fine ghost cells). It goes into the energies as the viscous work does,
  each phase taking its mass fraction of it.

Mass, momentum and total energy are conservative to round-off, across AMR levels as well: the
fluxes are refluxed (or, without subcycling, taken from the finer level at every stage), and in
`SIXEQS` the share of phase 2 is corrected with them. An absent
phase ($\alpha_k = 0$) gets no viscous energy and no heat, and a uniform velocity gives no
viscous terms.

In `SIXEQS_IE_NPHASE` the phase internal energies get the cell dissipation $\alpha_k\Phi_k$ and
$\rho E$ gets the conservative flux, so the sum of the phase internal energies and
$\rho E - \tfrac12\rho\lvert\vec u\rvert^2$ agree to truncation error only. With
`Physics.pressure_relaxation = 1` the phase internal energies are reset from $\rho E$ after
every stage.

### Material interfaces

Across an interface parallel to a face, the viscous traction and the heat flux are continuous
and the normal gradients of velocity and temperature jump, so the two half cells act in series.
The harmonic mean of the cell coefficients is the conductance of that series (Patankar 1980,
*Numerical Heat Transfer and Fluid Flow*, Sec. 4.2): with it a layered shear flow or a steady
heat flux through a sharp interface that lies on a cell face is exact, while the arithmetic mean
overestimates the face coefficient and gives a first-order error at the interface. In a uniform
mixture the two cells have the same coefficient and the harmonic mean is that coefficient exactly,
so mixture results do not change. A negative cell coefficient (from a volume fraction outside
$[0, 1]$) counts as 0 in the mean. The test cases `Couette2Layer-*` (a two-layer Couette flow) and
`HeatConduction2Mat-*` (conduction between two materials) check this against the exact
solutions ([Two-layer Couette flow](../verification.md#two-layer-couette-flow),
[Conduction between two materials](../verification.md#conduction-between-two-materials)).

An interface inside a cell leaves a mixed cell with the mixture coefficient, and the first-order
error comes back with the size of the arithmetic rule: in `Couette2Layer-6Eq` with the interface in
the middle of a cell (variant `mid-cell`) the wall stress is off by 2.1, 1.0 and 0.5 % on 32, 64
and 128 cells. In a diffuse-interface computation the interface is almost never on a face; with
the one-cell tanh profile of `Couette2Layer-*` (variant `diffuse`) the harmonic mean lowers the
wall-stress error by about 7 % of itself (9.3 % instead of 10.0 % on 32 cells).

A diffuse interface is a mixture: its cells have the mixture coefficients
$\sum_k\alpha_k\mu_k$ and $\sum_k\alpha_k\kappa_k$ of the model, which across the interface
exceed the series value. A layered flow across a diffuse interface therefore converges at first
order, with an error that grows with the width of the interface in cells and with the ratio of
the coefficients of the materials.

**Known limitations.**

- **No heat exchange between phases without the pressure-temperature relaxation.** A phase
  conducts only through itself, and heat passes from one phase to another only through the
  relaxation. With the pressure relaxation alone, and always in `SIXEQS_IE_NPHASE`, which has no
  temperature relaxation, two materials at different temperatures that meet at an interface do
  not conduct heat to each other: heat crosses only through the mixed cells of the interface,
  less and less as the grid is refined, and not at all through a sharp interface. The remedy
  would be the N-phase pressure-temperature relaxation of Huang (2026, Sec. 3.1.2) in
  `SIXEQS_IE_NPHASE`, after which the mixture heat flux above applies; it is not implemented. In
  `HeatConduction2Mat-6Eq-N` the heat through the interface is 52, 28 and 14 % of the exact value
  on 128, 256 and 512 cells. A run with a nonzero conductivity and no temperature relaxation
  prints a warning at startup.
- The five-equation models average the mixture viscosity and conductivity arithmetically to the
  faces, with the first-order interface error above.

The time step is limited by `run.vnn` ([Time stepping](inputs.md#time-stepping)) with the
diffusivity of each face: the larger of $(\tfrac43\mu + \mu_B)/\rho$, with the face viscosities
above and $\rho$ the smaller of the two cell densities, and the thermal diffusivity. For the
phase heat fluxes that is $(w_k/\min\alpha_k)\,\kappa_k/(\rho_kc_{v,k})$ over the phases that
conduct through the face, with $\rho_k$ the smaller of the two phase densities (the factor
$w_k/\min\alpha_k$ lies between 1 and 2); for the mixture heat flux it is $\kappa_f$ over the
smaller $\sum_k\alpha_k\rho_kc_{v,k}$ of the two cells.

The walls are those of the five-equation models: a no-slip wall (`5`) has zero velocity, and
every wall is adiabatic, for each phase.

## Comparing the models

With a single fluid every model reduces to the Euler equations. The models differ only in mixture
cells. The six-equation models reach pressure equilibrium through relaxation instead of the
$K_k\,\nabla\cdot\vec{u}$ term, which makes them more robust for strong shocks and expansions
through mixtures. The five-equation model is cheaper. The Sod shock tube computed with all four
models is on [Verification and validation](../verification.md#sod-shock-tube).

## Phase-Field interface regularization

Numerical diffusion thickens a moving interface over time. `PhaseField.phase_field = 1` adds
terms that hold it at a fixed width $\eta$. Volume moves between the phases with fluxes
$\vec J_k$, $\sum_k \vec J_k = 0$, and each phase carries its own mass, momentum and energy
with it. The terms added to the right-hand sides of the model are

$$
\begin{aligned}
\frac{\partial \alpha_k}{\partial t} + \dots &= \nabla\cdot\vec{J}_k\\
\frac{\partial (\alpha_k\rho_k)}{\partial t} + \dots &= \nabla\cdot(\rho_k\vec{J}_k)\\
\frac{\partial (\rho \vec{u})}{\partial t} + \dots &= \nabla\cdot(\vec{m} \otimes \vec{u}), \qquad \vec m = \sum_k \rho_k \vec J_k\\
\frac{\partial (\rho E)}{\partial t} + \dots &= \nabla\cdot\Big(\tfrac12\lvert\vec u\rvert^2\,\vec m + \sum_k \rho_k e_k \vec{J}_k \Big)
\end{aligned}
$$

In the six-equation models each phase energy gets its own term,
$\nabla\cdot(\rho_k E_k\vec J_k)$ in `SIXEQS` or $\nabla\cdot(\rho_k e_k \vec J_k)$ in
`SIXEQS_IE_NPHASE`.

The terms are in divergence form, so mass, momentum and energy are conserved, and consistent, so
uniform velocity, pressure and temperature stay uniform across the interface (Huang and Johnsen
2023 and 2024 for the two- and N-phase five-equation models, Huang 2026 for the six-equation
models). Their effect on advected interfaces, against MUSCL-MC alone and THINC, is shown under
[Interface advection](../verification.md#interface-advection).

## Phase-Field mechanisms

`PhaseField.Mechanism` chooses the volume flux. The two-phase models `FIVEEQS` and `SIXEQS`
carry only $\alpha_1$, apply the mechanism to $\phi = \alpha_1$ and set
$\vec J_2 = -\vec J_1$. The N-phase models apply it to every $\alpha_k$, with $\phi = \alpha_k$,
and then combine the fluxes so that $\sum_k \vec J_k = 0$, which with two phases again gives
$\vec J_2 = -\vec J_1$, since every flux below changes sign when $\phi$ is replaced by $1-\phi$. The interface normal is $\vec n = \nabla\phi/\lvert\nabla\phi\rvert$ and the curvature
is $\kappa = -\nabla\cdot\vec n$. `ACDI` computes the same normal from $\psi$, defined in its row, as
$\vec n = \nabla\psi/\lvert\nabla\psi\rvert$; $\psi$ varies almost linearly across the
interface, so its differences give a more accurate discrete normal. The CAC and Cahn–Hilliard mechanisms use the double-well potential
$g(\phi) = \phi^2(1-\phi)^2$, and the CAC mechanisms also use the weight $W(\phi) = \phi(1-\phi)$.
`CDI` and `ACDI` use neither.

| `PhaseField.Mechanism` | Volume flux | Notes |
|---|---|---|
| `CAC-Adv` (default) | $\nabla\cdot\vec J = M\left[\nabla^2\phi - \dfrac{g'(\phi)}{\eta^2} + \lvert\nabla\phi\rvert\,\kappa\right] + L\,W(\phi)$ | conservative Allen–Cahn, $L$ makes $\int\nabla\cdot\vec J\,d\Omega = 0$ |
| `CAC-VanDerWaals` | the same without $\lvert\nabla\phi\rvert\,\kappa$ | the interface also moves by its curvature |
| `CDI` | $\vec J = M\left[\nabla\phi - \dfrac{\sqrt2}{\eta}\,\phi(1-\phi)\,\vec n\right]$ | conservative diffuse interface |
| `ACDI` | $\vec J = M\left[\nabla\phi - \dfrac{1}{2\sqrt2\,\eta}\,\mathrm{sech}^2\!\left(\dfrac{\psi}{\sqrt2\,\eta}\right)\vec n\right]$ | accurate CDI, with $\psi = \sqrt2\,\eta\,\mathrm{artanh}(2\phi-1)$ |
| `CH-Adv`, `CH-VanDerWaals` | $\vec J = M\nabla\mu$, $\mu = \dfrac{g'(\phi)}{\eta^2} - \nabla^2\phi$ $\left(- \lvert\nabla\phi\rvert\,\kappa\right)$ | Cahn–Hilliard; the bracketed term is in `CH-Adv` only |

The variants discretize the same fluxes differently.

- `CAC-Adv-NormalFace`, `CDI-NormalFace`, `ACDI-NormalFace` and `CH-Adv-NormalFace` compute
  $\vec n$ on the cell faces, from $\psi$ in `ACDI-NormalFace` and from $\phi$ in the others. The
  other mechanisms compute $\vec n$ at cell centers and interpolate it to the faces.
- `CDI-CompressionInterpolation` forms the compression term $\phi(1-\phi)\,\vec n$ at cell
  centers and interpolates it to the faces, where `CDI` interpolates $\phi$ and $\vec n$ to
  the faces first.

All mechanisms share the planar equilibrium profile
$\phi = \tfrac12\left[1 + \tanh\left(\psi/\sqrt2\,\eta\right)\right]$, with $\psi$ as in the
`ACDI` row.

- $\eta$ is `PhaseField.Eta_Multiplier` times the cell size at `amr.max_level`, the largest over
  the directions, unless `PhaseField.Eta` gives it directly. With AMR, keep `Eta_Multiplier`
  at 0.5 or below for the implicit `CAC-Adv` mechanism (see [Phase-Field](inputs.md#phase-field)).
- $M$ is set from $\eta$ and $\max\lvert\vec u\rvert$: $M = \tfrac12\,\eta\max\lvert\vec u\rvert$,
  and $M = 3\,\eta^3\max\lvert\vec u\rvert$ for the Cahn–Hilliard mechanisms. In the CAC and
  Cahn–Hilliard mechanisms $M$ is also multiplied by `PhaseField.Lambda`, which is 1 by default.
- The terms are applied in a separate step after each hyperbolic time step, implicitly by
  default through AMReX's multigrid solvers. `PhaseField.ID_ExplicitRC = 1` treats the
  flux-based mechanisms (CDI, ACDI and CH) explicitly. The CAC mechanisms are implicit only.

The parameters of this step are listed under [Phase-Field](inputs.md#phase-field). With the
bound-preserving limiter of the five-equation models, `FiniteVolume.ID_Bound = 1`, the volume
fractions also stay bounded. In the six-equation models the `WENO3B` and `WENO5B` schemes bound
the volume fractions and partial densities.
The papers for the limiter, for adaptive mesh refinement and for the WENO schemes are listed in
[Citing COMPAS](../about/citing.md).
