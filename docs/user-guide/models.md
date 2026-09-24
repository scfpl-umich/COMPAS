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
\rho_k e_k = A_k P_k + B_k + \rho_k D_k ,
$$

which covers the three equations of state in COMPAS. The run-time parameter `prob.EOS.EOS`
selects one for all phases.

| `prob.EOS.EOS` | $A_k$ | $B_k$ | $D_k$ |
|---|---|---|---|
| `0`, ideal gas | $1/(\gamma_k-1)$ | $0$ | $0$ |
| `1`, stiffened gas | $1/(\gamma_k-1)$ | $\gamma_k P_{\infty,k}/(\gamma_k-1)$ | $0$ |
| `2`, Mie–Grüneisen | $1/\Gamma_k$ | $-p_{\mathrm{ref},k}/\Gamma_k$ | $e_{\mathrm{ref},k}$ |

The constants $\gamma_k$, $P_{\infty,k}$, $\Gamma_k$, $p_{\mathrm{ref},k}$,
$e_{\mathrm{ref},k}$ and $c_{v,k}$ are the inputs `gamma`, `pinf`, `gGamma`, `pref`, `eref`
and `cv`, listed under [Equation of state and materials](inputs.md#equation-of-state-and-materials).

Where a temperature is needed, it comes from

$$
\rho_k e_k = \rho_k c_{v,k} T_k + \frac{B_k}{1+A_k} + \rho_k D_k .
$$

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
P = \frac{\rho e - \sum_k \alpha_k \left(B_k + \rho_k D_k\right)}{\sum_k \alpha_k A_k} .
$$

`FIVEEQS` has two phases. `FIVEEQS_NPHASE` has $N$ phases, set with `-DNPHASE=`$N$ in the
`GNUmakefile`.

`Physics.source_term` sets the volume-fraction source.

- `0` (default): $K_k = 0$, so $\alpha_k$ is only advected (Allaire et al. 2002).
- `1`: $K_k = \dfrac{\rho c^2}{\rho_k c_k^2} - 1$, with
  $\dfrac{1}{\rho c^2} = \sum_j \dfrac{\alpha_j}{\rho_j c_j^2}$ (Kapila et al. 2001).

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
  $\alpha_k\rho_ke_k = A_k\,\alpha_kP_k + \alpha_kB_k + \alpha_k\rho_kD_k$.
- The two energy equations add up to the conservative equation for $\rho E$, so total energy
  is conserved.
- The relaxation $\mathcal{R} = \mu_P(P_1 - P_2)$ is instantaneous, $\mu_P \to \infty$, and
  is applied after every Runge–Kutta stage.
- `Physics.pressure_relaxation = 1` relaxes to $P_1 = P_2$, with $P_I$ the relaxed pressure.
- `Physics.pressure_temperature_relaxation = 1` relaxes to $P_1 = P_2$ and $T_1 = T_2$.
- With neither, the phases keep their own pressures. Both are off by default.
- The model is solved without viscous and conductive fluxes and is built with
  `-DDIFFUSION=false`.

**Reference.** Z. Huang, A consistent and conservative Phase-Field method for compressible
multiphase flows with the six-equation model, arXiv:2609.18085 (2026),
[doi:10.48550/arXiv.2609.18085](https://doi.org/10.48550/arXiv.2609.18085).

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
- `Physics.pressure_relaxation = 1`: after every stage a Newton iteration finds the common
  pressure $P$ from

  $$
  \sum_k \frac{\alpha_k\,(P_k - P)}{\rho_k c_k^2(P)} = 0 ,
  $$

  the volume fractions follow with $P_I = P$, and the phase internal energies are reset from
  $\rho e = \rho E - \tfrac12\rho\lvert\vec u\rvert^2$. It is off by default.
- $N$ is set with `-DNPHASE=`$N$ in the `GNUmakefile`.
- As for `SIXEQS`, the model is built with `-DDIFFUSION=false`.

**Reference.** The same paper as `SIXEQS`: Z. Huang, arXiv:2609.18085 (2026),
[doi:10.48550/arXiv.2609.18085](https://doi.org/10.48550/arXiv.2609.18085).

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
