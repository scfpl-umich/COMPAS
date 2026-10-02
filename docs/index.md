```{raw} html
<pre class="compas-logo" role="img" aria-label="COMPAS">
   __________  __  _______  ___   _____
  / ____/ __ \/  |/  / __ \/   | / ___/
 / /   / / / / /|_/ / /_/ / /| | \__ \
/ /___/ /_/ / /  / / ____/ ___ |___/ /
\____/\____/_/  /_/_/   /_/  |_/____/</pre>
```

# COMPAS

**COmpressible Multiphase multiPhysics Adaptive Solver**, from the Scientific Computing and Flow
Physics Laboratory at the University of Michigan.

COMPAS simulates compressible multiphase flows with shocks, interfaces and large density ratios
on block-structured adaptive meshes. It is built on [AMReX](https://amrex-codes.github.io/amrex/),
which provides the grids, the parallelism and the I/O, and the same source runs on a laptop, a CPU
cluster or GPUs.

Its main method is consistent and conservative Phase-Field interface regularization for
diffuse-interface compressible models: the interface stays thin and bounded while mass, momentum
and energy are conserved, for two or N phases and with adaptive mesh refinement. The methods are
described in the papers listed under [Citing COMPAS](about/citing.md).

```{raw} html
<figure>
<video controls preload="metadata" playsinline muted loop width="720" style="max-width:100%;height:auto"
       poster="_static/cylrm3_poster.png">
  <source src="_static/cylrm3.mp4" type="video/mp4">
  <a href="_static/cylrm3.mp4">Download the video (MP4)</a>
</video>
</figure>
```

*A Mach 2 cylindrical shock converges on an SF$_6$ shell (teal) around a helium core (dark blue)
in air (sand), with THINC keeping the interfaces between the three materials sharp. Lower half:
pressure-gradient schlieren of the shocks and waves, the two interfaces (black), and the outlines
of the regions covered by AMR levels 1 to 4, which follow the interfaces and the shocks. See
[`exec/CylindricalRM-5Eq-N`](examples.md#cylindricalrm-5eq-n).*

## What it can do

- **Models.** Two-phase and N-phase five-equation models, a two-phase six-equation model with
  pressure or pressure-temperature relaxation, and an N-phase six-equation model with pressure
  relaxation. The model is chosen at compile time.
- **Thermodynamics.** Ideal gas, stiffened gas, Mie-Grüneisen and Noble-Abel stiffened-gas equations
  of state. Viscous stresses and heat conduction in every model, per phase in the six-equation
  models. User-defined source terms.
- **Interfaces.** Phase-Field regularization with the Conservative Allen-Cahn, Cahn-Hilliard,
  Conservative Diffuse Interface and Accurate Conservative Diffuse Interface mechanisms, treated
  explicitly or implicitly. THINC interface reconstruction.
- **Numerics.** Finite volumes with Godunov, MUSCL, WENO3 and WENO5 reconstruction and bounded
  WENO variants, LLF, HLL and HLLC Riemann solvers, Forward Euler, TVD-RK2, TVD-RK3 and RK4
  time integration, and an optional bound-preserving limiter for the five-equation models.
- **Adaptive mesh refinement.** Refinement on gradients, thresholds or ranges of any variable,
  optional subcycling in time, and conservative refluxing.
- **Parallelism.** MPI, OpenMP, and GPUs through CUDA, HIP or SYCL.
- **Geometry and I/O.** 2D and 3D Cartesian domains, STL geometry for initial conditions, AMReX
  plotfiles readable by ParaView, VisIt and yt, and checkpoint and restart.

The methods are chosen at run time from a text inputs file, so comparing them needs no rebuild.
The options are listed in [Compile-time options](user-guide/compile-options.md) and
[Run-time options](user-guide/inputs.md).

## Where to start

- [Installation](getting-started/install.md): requirements, getting COMPAS and AMReX, building,
  and the Python environment for analysis.
- [First simulation](getting-started/quickstart.md): build and run a small case, look at the
  results, and switch interface methods.
- [Models and equations](user-guide/models.md): the equations each model solves.
- [Setting up a case](user-guide/cases.md): the files of a case and how to write your own.
- [Examples](examples.md): the cases provided in `exec/` and `exec/_Tests/`.
- [Testing](user-guide/testing.md): the test suite to run before opening a pull request.
- [Verification and validation](verification.md): comparisons with exact solutions and
  experiments.

COMPAS is free software under the GNU General Public License, version 3 or later. See
[License](about/license.md) and [Contributing](about/contributing.md).

```{toctree}
:hidden:
:caption: Getting started

getting-started/install
getting-started/quickstart
```

```{toctree}
:hidden:
:caption: User guide

user-guide/models
user-guide/cases
user-guide/compile-options
user-guide/inputs
user-guide/testing
```

```{toctree}
:hidden:
:caption: Examples and checks

examples
verification
```

```{toctree}
:hidden:
:caption: About

about/citing
about/license
about/contributing
```
