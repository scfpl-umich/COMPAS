# COMPAS

**COmpressible Multiphase multiPhysics Adaptive Solver**
Scientific Computing and Flow Physics Laboratory, University of Michigan

COMPAS simulates compressible multiphase flows with shocks, interfaces and large density
ratios on block-structured adaptive meshes. It is built on [AMReX](https://amrex-codes.github.io/amrex/)
and runs on laptops, CPU clusters and GPUs from the same source.

Its main method is **consistent and conservative Phase-Field interface regularization** for
diffuse-interface compressible models: the interface stays sharp and bounded while mass, momentum
and energy are conserved, for two or N phases and with adaptive mesh refinement. The methods are
described in the papers listed under [Citing COMPAS](#citing-compas).

<p align="center">
  <img src="https://raw.githubusercontent.com/scfpl-umich/COMPAS-media/main/figures/cylrm3_hero.png" width="420"
       alt="Three-material cylindrical Richtmyer-Meshkov implosion computed with COMPAS">
</p>

*`exec/CylindricalRM-5Eq-N` at t = 0.4: a converging shock drives SF6 jets out of a shell (teal)
around a helium core (dark blue) in air (sand), with THINC interfaces. Lower half: pressure-gradient
schlieren, the two interfaces (black), and the regions covered by AMR levels 1 to 4, which follow
the interfaces and the shocks. The scale bar is 0.1.*

## Features

**Multiphase models** (selected at compile time)

| `PHYSICS` | Model |
|---|---|
| `FIVEEQS` | Two-phase five-equation model, with or without the compressibility source term |
| `FIVEEQS_NPHASE` | N-phase five-equation model (`NPHASE` set at compile time) |
| `SIXEQS` | Two-phase six-equation model with pressure or pressure–temperature relaxation |
| `SIXEQS_IE_NPHASE` | N-phase six-equation model, internal-energy formulation, with pressure relaxation |

**Thermodynamics and physics**
- Equations of state: ideal gas, stiffened gas, Mie–Grüneisen (constant Γ with reference pressure and energy), and Noble-Abel stiffened gas
- Viscous stresses and heat conduction (all models; per phase in the six-equation models)
- User-defined source terms through a hook in each case's `ProblemICBC.H`

**Phase-Field interface regularization**
- Conservative Allen–Cahn (CAC), Cahn–Hilliard (CH), Conservative Diffuse Interface (CDI) and
  Accurate Conservative Diffuse Interface (ACDI) mechanisms, including face-normal variants
- Explicit or implicit treatment, the latter through AMReX's multigrid linear solvers
- Consistent with the multiphase model and conservative, with bounded volume fractions when
  the bound-preserving limiter of the five-equation models is enabled

**Numerics**
- Finite-volume discretization with Godunov, MUSCL (MC, minmod, superbee, van Leer), WENO3 and
  WENO5, and their bounded variants WENO3B and WENO5B, for every model
- THINC interface reconstruction, and a quadrature-based high-order WENO5 (`FIVEEQS`, 2D)
- Local Lax–Friedrichs, HLL and HLLC Riemann solvers
- Forward Euler, TVD-RK2, TVD-RK3 and RK4 time integration
- Optional bound-preserving limiter for the five-equation models (`FIVEEQS`, `FIVEEQS_NPHASE`)
- Finite-difference operators for gradients, interface normals, curvature and diffusive fluxes

**Adaptive mesh refinement and parallelism**
- Patch-based AMR with refinement on gradients, thresholds or ranges of any variable
- Temporal subcycling (optional), conservative refluxing, and piecewise-constant, linear or
  quartic coarse–fine interpolation
- MPI, OpenMP, and GPUs (CUDA, HIP, SYCL) through AMReX

**Geometry and I/O**
- 2D and 3D Cartesian domains
- STL surface reader for initializing complex geometries
- AMReX plotfiles, readable directly by ParaView, VisIt and yt
- Checkpoint and restart, and time-windowed output

## Documentation

The documentation is at <https://scfpl-umich.github.io/COMPAS/>: installation, a first
simulation, the models and their equations, setting up a case, every compile-time and run-time
option, the example cases, and verification and validation.

## Requirements

- A C++17 compiler (GCC or Clang) and GNU Make
- MPI for parallel runs (optional; set `USE_MPI = FALSE` for a serial build)
- [AMReX](https://github.com/AMReX-Codes/amrex), cloned next to COMPAS. COMPAS has been
  tested with AMReX 26.01.
- Python 3, only for the optional convergence and plotting scripts in `scripts/`

## Quickstart

Clone COMPAS and AMReX side by side:

```bash
mkdir compas-workspace && cd compas-workspace
git clone https://github.com/scfpl-umich/COMPAS.git
git clone --branch 26.01 https://github.com/AMReX-Codes/amrex.git
```

Build a two-dimensional advected interface with two levels of refinement, then run it on four
MPI ranks (about 15 seconds on a laptop):

```bash
cd COMPAS/exec/_Tests/Advection-5Eq
make -j4
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs amr.plot_int=50
```

The executable name records the dimension, the compiler, and whether MPI, OpenMP or a GPU
backend is enabled. A serial build (`USE_MPI = FALSE`) runs without `mpirun`.

Plotfiles are written to `plot/Advection5Eq/` (`amr.case_name` in the inputs file). Open them in
ParaView or VisIt as AMReX plotfiles, or with [yt](https://yt-project.org/); the `.visit` file in
the same folder loads the whole run as a time series in VisIt. `bash ../../clean.sh` removes the
output of previous runs.

The same problem with Phase-Field interface regularization is in `exec/_Tests/Advection-PF-5Eq`
(the default CAC-Adv mechanism, about a minute on four ranks), and with THINC in
`exec/_Tests/Advection-THINC-5Eq`. Both build and run the same way.

Any parameter in `prob/inputs` can be overridden on the command line, the quickest way to try a
different method:

```bash
./main2d.gnu.MPI.ex prob/inputs FiniteVolume.Scheme=WENO5 Physics.RiemannSolver=HLL
```

## Setting up a case

A case is a directory with a `GNUmakefile` and a `prob/` folder holding `Parm.H` (problem
parameters), `ProblemICBC.H` (the initial condition and optional hooks) and `inputs` (run-time
parameters). Start from the closest case in `exec/_Tests/` or `exec/` and copy it.

The documentation has the details:

- [Setting up a case](https://scfpl-umich.github.io/COMPAS/user-guide/cases.html): the files, the hooks and in-situ diagnostics
- [Compile-time options](https://scfpl-umich.github.io/COMPAS/user-guide/compile-options.html): the model, dimension and backends
- [Run-time options](https://scfpl-umich.github.io/COMPAS/user-guide/inputs.html): every inputs parameter and its default
- [Examples](https://scfpl-umich.github.io/COMPAS/examples.html): the cases provided in `exec/`

## Repository layout

```
COMPAS/
├── source/        driver (AMReX AmrCore) and model classes
├── include/       physics, operators and tools, header-only for GPU portability
│   ├── Physics_*          model definitions: state, fluxes, equations of state
│   ├── Operators_*        discretization: finite volume, finite difference,
│   │                      phase field, linear systems
│   └── Tools/             STL reader and utilities
├── exec/
│   ├── _Tests/    small cases run by the test suite
│   └── ...        larger application cases
└── scripts/       test runner, convergence and plotting tools
```

## Testing

```bash
bash scripts/test_cases.sh
```

The suite builds every case in `exec/_Tests/` and briefly runs it and the variants in its
`variants.txt`, checking each run for its exit status, for finite values in its log, and for the
use of every override it is given. It is a build-and-stability check; the comparisons with
reference solutions are under
[Verification and validation](https://scfpl-umich.github.io/COMPAS/verification.html), and the
runner's options under [Testing](https://scfpl-umich.github.io/COMPAS/user-guide/testing.html).

## Contributing

Contributions are welcome through pull requests against `main`, which a maintainer reviews and
approves before merging. [CONTRIBUTING.md](CONTRIBUTING.md) describes what a pull request needs,
including the test suite and what to do when a change affects the cases.

## Citing COMPAS

If you use COMPAS, please cite the papers that describe the methods you used.

- **COMPAS.** W. J. White, Z. Huang, and E. Johnsen, COMPAS: A GPU-accelerated, adaptive-mesh finite volume framework for consistent and conservative Phase-Field modeling of compressible, N-material, interfacial flows, in preparation (2026).
- **COMPAS, PhD dissertation.** W. J. White, Discontinuous Galerkin Methods and Phase-Field Models for Numerical Simulations of Compressible Interfacial Flows with Shocks and Vorticity, PhD dissertation, University of Michigan (2025), [doi:10.7302/28171](https://doi.org/10.7302/28171). One of its chapters presents COMPAS and its verification and validation, and is the basis of the COMPAS paper above.

- **Two-phase five-equation model with Phase-Field.** Z. Huang and E. Johnsen, A consistent and conservative Phase-Field method for compressible multiphase flows with shocks, *Journal of Computational Physics* 488 (2023) 112195, [doi:10.1016/j.jcp.2023.112195](https://doi.org/10.1016/j.jcp.2023.112195).
- **N-phase five-equation model with Phase-Field.** Z. Huang and E. Johnsen, A consistent and conservative Phase-Field method for compressible N-phase flows: Consistent limiter and multiphase reduction-consistent formulation, *Journal of Computational Physics* 501 (2024) 112801, [doi:10.1016/j.jcp.2024.112801](https://doi.org/10.1016/j.jcp.2024.112801).
- **Bound preservation.** Z. Huang and E. Johnsen, Bound preservation for the consistent and conservative Phase-Field method for compressible single-, two-, and N-phase flows, *Journal of Computational Physics* 526 (2025) 113783, [doi:10.1016/j.jcp.2025.113783](https://doi.org/10.1016/j.jcp.2025.113783).
- **Adaptive mesh refinement.** Z. Huang, W. J. White, and E. Johnsen, Consistent and conservative Phase-Field method for compressible two- and N-phase flows with adaptive mesh refinement, *Journal of Computational Physics* 548 (2026) 114569, [doi:10.1016/j.jcp.2025.114569](https://doi.org/10.1016/j.jcp.2025.114569).
- **Consistent and bound-preserving WENO.** Z. Huang, Consistent and bound-preserving finite-volume WENO scheme for compressible two-/N-phase flows with Phase-Field mechanism, arXiv:2608.00746 (2026), [doi:10.48550/arXiv.2608.00746](https://doi.org/10.48550/arXiv.2608.00746).
- **Six-equation model with Phase-Field.** Z. Huang, A consistent and conservative Phase-Field method for compressible multiphase flows with the six-equation model, arXiv:2609.18085 (2026), [doi:10.48550/arXiv.2609.18085](https://doi.org/10.48550/arXiv.2609.18085).

## Acknowledgments

This work was supported in part by Los Alamos National Laboratory under the project “Algorithm/Software/Hardware Co-design for High Energy Density applications” at the University of Michigan.

## License

COMPAS is free software, distributed under the GNU General Public License, version 3 or later
([LICENSE](LICENSE)). Copyright (C) 2026 The Regents of the University of Michigan. The contributors
are listed in [AUTHORS](AUTHORS).

[COPYRIGHT](COPYRIGHT) adds one permission to the license: you may convey COMPAS linked with the
vendor GPU toolkits and drivers (NVIDIA, AMD), the Intel oneAPI and HPE Cray programming
environments, any MPI implementation, and compiler runtime libraries, including proprietary ones.
Parts of COMPAS are derived from AMReX, whose BSD-3-Clause license and notice are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Running COMPAS, modifying it for your own use and publishing results obtained with it place no
obligation on you under the license. Its conditions apply when you convey COMPAS or a modified
version to others, in source or binary form, including executables and container images. Citing
the papers above is requested but is not a condition of the license.
