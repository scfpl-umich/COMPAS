# Setting up a case

A case is a directory that holds everything specific to one problem: the compile-time choices,
the problem parameters, the initial condition and any user hooks, and the run-time inputs. The
solver itself lives in `source/` and `include/` and is shared by every case.

## The files in a case

```
MyCase/
├── GNUmakefile
└── prob/
    ├── Make.package
    ├── Parm.H
    ├── ProblemICBC.H
    └── inputs
```

| File | Contents |
|---|---|
| `GNUmakefile` | Compile-time choices: the model, the dimension and the parallel backends |
| `prob/Make.package` | The problem headers for the build. Rarely changed |
| `prob/Parm.H` | `struct Parm`, the problem's parameters and the EOS constants, read from the `prob.*` entries of the inputs file |
| `prob/ProblemICBC.H` | The initial condition, and optional hooks for source terms, boundary conditions, refinement and output |
| `prob/inputs` | Everything chosen at run time: domain, grid, methods, output and problem values |

## Starting a new case

Copy the closest existing case in `exec/` or `exec/_Tests/` and edit it. The
[Examples](../examples.md) page lists them. For example,

```bash
cd COMPAS/exec
cp -r ShockBubble-AirHelium-5Eq MyCase
cd MyCase
make -j4
```

A copy in `exec/` or `exec/_Tests/` builds as is, because its `GNUmakefile` ends with
`include ../Make.exec`, which finds the build rules. For a case kept outside the repository,
change that line to

```make
include $(COMPAS_HOME)/exec/Make.exec
```

and build with

```bash
make -j4 COMPAS_HOME=/path/to/COMPAS
```

AMReX is expected next to COMPAS. `AMREX_HOME=/path/to/amrex` points the build at another
checkout.

Change `amr.case_name` in the new `prob/inputs`, so that the copy writes its output to its own
folder. `bash ../clean.sh` (from a case in `exec/`) or `bash ../../clean.sh` (from a case in
`exec/_Tests/`) removes the plotfiles and checkpoints of previous runs.

## GNUmakefile

The `GNUmakefile` sets the model (`-DPHYSICS=`), the number of phases for the N-phase models
(`-DNPHASE=`), the dimension (`DIM`) and the backends (`USE_MPI`, `USE_OMP`, `USE_CUDA`,
`USE_HIP`, `USE_SYCL`). Every option is described in [Compile-time options](compile-options.md).

## Make.package

`prob/Make.package` lists the headers of the case:

```make
CEXE_headers += ProblemICBC.H
CEXE_HEADERS += Parm.H
```

Add a line for any new header the case includes, as `exec/GuderleyImplosion-5Eq` does for
`prob/Guderley.H`.

## Parm.H

`Parm.H` defines `struct Parm`, which holds the parameters of the problem. COMPAS fills it on the
host at startup and copies it to the GPU, and every kernel and hook receives it as `parm`.

`Parm` contains the parameter structures of the physics, finite-volume, Phase-Field and linear
solver modules, and the EOS and transport constants of each material. The member function
`DynamicInit()` reads the constants from `prob.EOS.*` and then calls `Initialize()`, which reads
the `Physics.*`, `FiniteVolume.*` and `PhaseField.*` entries. The two-phase models name the
constants per material (`prob.EOS.gamma_1`, `prob.EOS.pinf_2`, ...), and the N-phase models read
one list per property (`prob.EOS.gamma = 50.0 1.66667 1.4`). The default EOS of a case is the
initial value of `EOS` in its `Parm.H`.

To add a problem parameter, add a member with a default value and read it from `prob.*` at the
end of `DynamicInit()`. From `exec/_Tests/RichtmyerMeshkov-5Eq/prob/Parm.H`:

```cpp
struct Parm
{
    // ... physics, EOS and solver members ...

    amrex::Real rho1 = 1.0;    // light-gas (material 1) density ahead of the shock
    amrex::Real Mach = 1.5;    // shock Mach number in material 1

    void Initialize();
    #define DYNAMIC_INIT
    void DynamicInit();
};

AMREX_FORCE_INLINE
void Parm::DynamicInit ()
{
    amrex::ParmParse pp("prob.EOS");
    // ... EOS constants ...
    Initialize();

    amrex::ParmParse ppp("prob");
    ppp.query("rho1", rho1);
    ppp.query("Mach", Mach);
}
```

The member is then `prob.Mach` in the inputs file and `parm.Mach` in the kernels (`parm->Mach`
in the boundary functions, which receive a pointer).

`Parm` is copied to the GPU byte for byte, so its members must be plain values and fixed-size
arrays. An array whose size is known only at run time goes in device memory, with a pointer to
it in `Parm`. `exec/_Tests/RichtmyerMeshkov-Multimode-5Eq/prob/Parm.H` does this for its list of
interface modes.

## ProblemICBC.H

### Initial condition

Every case defines `initial_condition`:

```cpp
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
```

It returns the conservative state at the point $(x, y, z)$. `dx` is the cell size of the level
being filled and `dx_FinestLevel` that of level `amr.max_level`, even before that level exists,
which cases use to set the width of a diffuse interface. The cases set `sharp = false`, and the
core does not read it.

By default the function is evaluated at 4 Gauss-Legendre points in each direction of every cell
and the results are averaged, so the initial data are cell averages. `run.QuadratureIC = 0`
evaluates it once, at the cell center.

The state is easiest to build from primitive variables with `State2Primitive` and
`Primitive2Conservative`, so that the same code works with any EOS. Their arguments depend on the
model (the six-equation models take a pressure for each phase), so copy the calls from a case of
the same model.

### Optional hooks

Each hook is enabled by defining a macro in `ProblemICBC.H` just before the function; without it,
the core uses its default behavior.

| Macro | Function | Called | Example |
|---|---|---|---|
| `USER_SOURCE_TERM` | `user_source_term` | In every cell at every stage of the time step | `exec/_Tests/RayleighTaylor-5Eq` |
| `USER_BOUNDARY_FUNC_ROBIN` | `BoundaryFunction_Robin` | Ghost cells of faces with boundary code 1 | `exec/_Tests/Jet-Inflow-5Eq` |
| `USER_BOUNDARY_FUNC_DIRICHLET`, `USER_BOUNDARY_FUNC_NEUMANN` | `BoundaryFunction_User` | Ghost cells of faces with boundary code 1 | |
| `USER_BOUNDARY_FUNC` | `BoundaryFill_User` | Ghost cells outside the domain | `exec/GuderleyImplosion-5Eq` |
| `SPATIO_TEMPORAL_TAG_USER_FUNC_` | `spatio_temporal_tag` | In every cell of each level being tagged, when the grids are rebuilt | `exec/_Tests/NonsphericalCollapse-6Eq` |
| `USER_OUTPUT_FUNC` | `UserOutputFunction` | See [In-situ diagnostics](#in-situ-diagnostics) | `exec/_Tests/RichtmyerMeshkov-5Eq` |
| `CONVERGENCE` | `exact_solution` | At the end of the run, for error norms | `exec/_Tests/IsentropicVortex-5Eq` |

`include/BoundaryCondition.H` holds templates of the three boundary functions with every argument
in place.

#### Source term

```cpp
#define USER_SOURCE_TERM
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void user_source_term (
    amrex::Array<amrex::Real, NSTATE>& Residual,
    amrex::Array4<const amrex::Real> const& statein,
    int i, int j, int k,
    amrex::Real x, amrex::Real y, amrex::Real z,
    amrex::Real t,
    amrex::Real dt,
    AMREX_D_DECL(amrex::Real dx,
                 amrex::Real dy,
                 amrex::Real dz),
    Parm const& parm)
```

`Residual` is added to the right-hand side of the conservative equations in cell `(i,j,k)`, whose
state is `statein(i,j,k,n)`. The array is zeroed before the call, so a component the function
leaves unset adds nothing. `exec/_Tests/RayleighTaylor-5Eq` adds gravity to the momentum and energy
equations this way.

#### Boundary conditions

Boundary code 1 (`run.lo_bc`, `run.hi_bc`) is an inflow whose ghost cells are filled by a user
function, so code 1 is used together with one of the boundary macros below.

With `USER_BOUNDARY_FUNC_ROBIN`, `BoundaryFunction_Robin` is called for each ghost cell
`(i,j,k)` and each state component `n`, and sets three coefficients of the condition

$$
a\, q_b + b\, \frac{\partial q}{\partial n}\Big|_b = f ,
$$

where $q_b$ is the value of component `n` on the boundary face and $\partial q/\partial n$
its outward normal derivative there:

```cpp
#define USER_BOUNDARY_FUNC_ROBIN
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_Robin(int i, int j, int k, int n,
                            Real& BC_Dirichlet,     // a
                            Real& BC_Neumann,       // b
                            Real& BC_Function,      // f
                            Real const& xBC,
                            Real const& yBC,
                            Real const& zBC,
                            Real const& time,
                            Array4<Real const> const& q,
               AMREX_D_DECL(int ilo, int jlo, int klo),
               AMREX_D_DECL(int ihi, int jhi, int khi),
               AMREX_D_DECL(Real const& dx,
                            Real const& dy,
                            Real const& dz),
                            Parm const* parm)
```

`(xBC, yBC, zBC)` is the point on the boundary face, and comparing `i` with `ilo` and `ihi` (and
likewise for `j` and `k`) tells which face the ghost cell is on. `exec/_Tests/Jet-Inflow-5Eq`
imposes a jet state through part of the $x = 0$ face and a slip wall on the rest.

`USER_BOUNDARY_FUNC_DIRICHLET` and `USER_BOUNDARY_FUNC_NEUMANN` are the two special cases. The
case defines `BoundaryFunction_User`, with the same arguments except that only `BC_Function` is
set, and it is the boundary value or the normal derivative.

With `USER_BOUNDARY_FUNC`, the case fills the ghost cells itself:

```cpp
#define USER_BOUNDARY_FUNC
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFill_User(const IntVect& iv,
                       Array4<Real> const& q,
                       const int dcomp,
                       const int numcomp,
                       GeometryData const& geom,
                       const Real time,
                       const BCRec* bcr,
                       const int bcomp,
                       const int orig_comp,
                       Parm const* parm)
```

It is called for ghost cells on every face, so it checks `bcr[bcomp]` for `BCType::ext_dir`
(code 1) before acting. It is also called for fields other than the state, such as the
Phase-Field buffers. `exec/GuderleyImplosion-5Eq` tells these apart with `dcomp`, `orig_comp` and
`numcomp`, and fills the state with the exact similarity solution at the current time.

#### Refinement in space and time

`spatio_temporal_tag` tags cells for refinement by position, level and time, in addition to the
`run.refine.*` criteria of the inputs file. `SPATIO_TEMPORAL_TAG_USER_FUNC_` replaces the empty
default function and makes the tagging routine call it. `exec/_Tests/NonsphericalCollapse-6Eq`
also defines `SPACE_TIME_REF`, which is not required.

```cpp
#define SPATIO_TEMPORAL_TAG_USER_FUNC_
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void spatio_temporal_tag(int& TagVal,
                         int lev,
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z,
                         amrex::Real t,
                         Parm const& parm)
{
    if (x < 1.25 && y < 1.25) { TagVal = TagBox::SET; }
}
```

`(x, y, z)` is the cell center. Leave `TagVal` unchanged to leave the cell to the other criteria.

#### Exact solution

A case that defines `CONVERGENCE` as `true` and provides `exact_solution` gets the error norms of
every conservative variable against it at the end of the run. `scripts/convergence.py` uses them;
see [Testing](testing.md#convergence-script).

```cpp
#undef  CONVERGENCE          // Macros.H defines it as false by default
#define CONVERGENCE true
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void exact_solution(amrex::Array<amrex::Real, NSTATE>& Conservative,
                    amrex::Real x, amrex::Real y, amrex::Real z,
                    amrex::Real t,
                    AMREX_D_DECL(amrex::Real dx, amrex::Real dy, amrex::Real dz),
                    AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                 amrex::Real dy_FinestLevel,
                                 amrex::Real dz_FinestLevel),
                    Parm const& parm)
```

## inputs

`prob/inputs` is a list of `name = value` lines grouped by prefix, and the header at its top
describes the case. Every entry can be overridden on the command line. The entries under `prob.`
are the ones read by the case's `Parm.H`. The rest are described in
[Run-time options](inputs.md).

## In-situ diagnostics

A quantity needed at every step, such as a shock position, a mixing width or a peak pressure, is
better computed in memory during the run than from plotfiles, which are written far less often. A
case does this in `UserOutputFunction(T& o)`, where `o` is the solver object, and defines
`USER_OUTPUT_FUNC` to enable it. COMPAS calls the function once before the first step, and then

- every N coarse steps with `run.user_output_int = N`, or
- with each plotfile with `run.user_output_int = 0`, the default.

It is also called once before a run stops on a NaN in the step sum. Every rank calls it, and the
helpers in `include/Tools/UserOutput.H` (namespace `UserOutput`) do the shared work, including
the MPI reductions and the file output.

### Measure

```cpp
Values<N> Measure(o, {op_1, ..., op_N}, ngrow, f);
Values<N> Measure<N>(o, op, ngrow, f);
Values<N> Measure(o, {op_1, ..., op_N}, ngrow, f, FinestLevel);
```

`Measure` evaluates the GPU lambda `f` in each cell of the composite grid, the cells of each
level that no finer level covers. `f` returns `N` values, and output `n` is combined over the
cells and the MPI ranks with `op_n`, which is `Sum`, `Min` or `Max`. The second form applies one
operation to every output. With `FinestLevel` as the last argument, `f` is evaluated on the
cells of the finest level only, which may cover only part of the domain. The result is a
`Values<N>`, an `amrex::GpuArray<Real, N>`, and is the same on every rank.

The lambda has the form

```cpp
[=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<N> { ... }
```

| Argument | Contents |
|---|---|
| `c.i`, `c.j`, `c.k` | the index of the cell on its level |
| `c.x`, `c.dx` | its center and size, one value per direction |
| `c.dV` | its volume, an area in 2D |
| `c.domain` | the index box of the domain on its level |
| `c.finest` | `true` on the finest level |
| `U` | the conservative state, `U(i,j,k,n)` with `n` one of the `INDEX_*` constants of the model |
| `parm` | the case's `Parm` |

A `Sum` adds the values themselves, so an integral returns its integrand times `c.dV`. A cell
with nothing to contribute returns 0 for a `Sum`, a large value such as `1e30` for a `Min`, and a
large negative value such as `-1e30` for a `Max`.

`ngrow` is the number of ghost cells of `U`. With `ngrow = 0`, `f` may read only `U` at
`(c.i, c.j, c.k)`. A lambda that reads neighbors, for a gradient or a crossing, sets `ngrow` to
the distance it reaches, and `Measure` fills those cells from the same level, the coarser levels
and the boundary conditions.

### WriteRow and the accessors

`WriteRow(o, file, {{"name_1", value_1}, ...})` appends one line with the time, the coarse step
and the values to `file` in the plot directory of the case, for example
`./plot/Advection5Eq/advection.csv`. Only the I/O rank writes. At step 0 it starts a new file with
the header `t,step,name_1,...`, and a restarted run appends to the existing file.

| Accessor | Returns |
|---|---|
| `Step(o)` | the coarse step |
| `Time(o)` | the time |
| `HostParm(o)` | the case's `Parm` on the host |
| `FinestGeom(o)` | the `amrex::Geometry` of the finest level |
| `Directory(o)` | the plot directory, where `WriteRow` writes |

Anything the helpers do not cover can use the solver's members directly, such as
`o.dof_new[lev]`, `o.Geom(lev)` and `o.finestLevel()`.

### Example

`exec/_Tests/ShuOsher-5Eq` records the shock position, the largest $x$ at which the pressure
crosses a threshold between a cell and its $+x$ neighbor. The neighbor needs one ghost cell, and
the crossing is taken on the finest level:

```cpp
#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real P_thr = 1.5 * 1.0;

    // Finest level only, with one ghost cell for the crossing with the +x neighbor
    const auto v = Measure(o, {Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<1>
    {
        amrex::Array<Real, NSTATE> Uc, Ue;
        for (int n = 0; n < NSTATE; ++n) { Uc[n] = U(c.i,c.j,c.k,n); Ue[n] = U(c.i+1,c.j,c.k,n); }
        Real Pc, Pe;
        Conservative2Pressure(Pc, Uc, parm);
        Conservative2Pressure(Pe, Ue, parm);
        const bool cross = (Pc - P_thr) * (Pe - P_thr) <= 0.0 && Pc != Pe;
        return {cross ? c.x[0] + c.dx[0] * (P_thr - Pc) / (Pe - Pc) : -1.0e30};
    }, FinestLevel);
    WriteRow(o, "shock.csv", {{"x_shock", v[0]}});
}
```

Every case in `exec/` and `exec/_Tests/` except `IsentropicVortex-5Eq`, which reports error
norms instead ([Exact solution](#exact-solution)), defines a `UserOutputFunction` with these
helpers, with a comment above it that describes what it writes.
