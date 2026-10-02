# Compile-time options

Compile-time options are set in the case's `GNUmakefile` ([Setting up a case](cases.md)) or on
the `make` command line, for example `make -j4 DIM=3 USE_OMP=TRUE`. Changing one requires a
rebuild. Everything else is chosen at run time ([Run-time options](inputs.md)).

## Model

The model is fixed at compile time because the number of state variables differs between
models. Each option is a line `DEFINES += -D<NAME>=<value>` in the `GNUmakefile`. When a line is
absent, the default in `include/Macros.H` applies.

| Option | Values | Default | Meaning |
|---|---|---|---|
| `DEFINES += -DPHYSICS=` | `FIVEEQS`, `FIVEEQS_NPHASE`, `SIXEQS`, `SIXEQS_IE_NPHASE` | `FIVEEQS` | multiphase model ([Models and equations](models.md)) |
| `DEFINES += -DNPHASE=` | integer | `2` | number of phases for the N-phase models |
| `DEFINES += -DADVECTION=` | `true` / `false` | `true` | hyperbolic fluxes |
| `DEFINES += -DDIFFUSION=` | `true` / `false` | `false` | viscous and conductive fluxes |
| `DEFINES += -DNONCONSERVATIVE=` | `true` / `false` | `true` | non-conservative terms of the model |

`ADVECTION` and `NONCONSERVATIVE` are normally left `true`. `DIFFUSION=true` applies to every
model; the six-equation models give each phase its own stress and heat flux
([Viscous and conductive fluxes in the six-equation models](models.md#viscous-and-conductive-fluxes-in-the-six-equation-models)).
`NPHASE` is used by `FIVEEQS_NPHASE` and `SIXEQS_IE_NPHASE`. The two-phase models ignore
it.

The number of state variables in $d$ dimensions is

| `PHYSICS` | State variables | Count |
|---|---|---|
| `FIVEEQS` | $\alpha_1\rho_1,\ \alpha_2\rho_2,\ \rho\vec u,\ \rho E,\ \alpha_1$ | $4+d$ |
| `FIVEEQS_NPHASE` | $\alpha_k\rho_k\ (k=1..N),\ \rho\vec u,\ \rho E,\ \alpha_k\ (k=1..N)$ | $2N+d+1$ |
| `SIXEQS` | $\alpha_1\rho_1,\ \alpha_2\rho_2,\ \rho\vec u,\ \alpha_1\rho_1E_1,\ \alpha_2\rho_2E_2,\ \alpha_1$ | $5+d$ |
| `SIXEQS_IE_NPHASE` | $\alpha_k,\ \alpha_k\rho_k,\ \alpha_k\rho_ke_k\ (k=1..N),\ \rho\vec u,\ \rho E$ | $3N+d+1$ |

The names of these variables in plotfiles and refinement criteria are listed under
[Output variables](inputs.md#output-variables).

## Build

These are AMReX build options, set as `NAME = value` lines in the `GNUmakefile`.

| Option | Values | Meaning |
|---|---|---|
| `DIM` | `2`, `3` | spatial dimension |
| `USE_MPI`, `USE_OMP` | `TRUE` / `FALSE` | parallel backends |
| `USE_CUDA`, `USE_HIP`, `USE_SYCL` | `TRUE` / `FALSE` | GPU backend |
| `DEBUG` | `TRUE` / `FALSE` | debug build with assertions and array bounds checks, run with checkpoints off (`amr.chk_int = -1`) |

The case `GNUmakefile`s also set `COMP = gnu`, `PRECISION = DOUBLE` and `PROFILE = FALSE`, which
are passed to the AMReX build system unchanged. How the executable name records these options is
described under [The executable name](../getting-started/install.md#the-executable-name).

## Paths

`exec/Make.exec`, which every case includes, sets two paths.

| Variable | Default | Meaning |
|---|---|---|
| `COMPAS_HOME` | the repository that contains `exec/Make.exec` | COMPAS source tree |
| `AMREX_HOME` | `$(COMPAS_HOME)/../amrex`, a checkout next to COMPAS | AMReX source tree |

Both can be set in the environment or on the `make` line. A case kept outside the repository
includes `$(COMPAS_HOME)/exec/Make.exec` and is built with `make COMPAS_HOME=/path/to/COMPAS`
([Starting a new case](cases.md#starting-a-new-case)).

## Macros in the case files

Some options are macros defined in the case's `prob/ProblemICBC.H` instead of the
`GNUmakefile`. They enable the optional hooks for source terms, boundary conditions,
refinement and output, described under [Optional hooks](cases.md#optional-hooks), and the
error norms of a convergence study (`CONVERGENCE`), described under
[Convergence studies](inputs.md#convergence-studies).
