# Installation

COMPAS is not installed as a library. Each case is its own small program: `make` in a case
directory compiles COMPAS, AMReX and the problem files together into one executable. Installing
COMPAS therefore means getting the two source trees and a compiler.

## Requirements

- A C++17 compiler (GCC or Clang) and GNU Make.
- MPI for parallel runs, for example MPICH or Open MPI. A serial build needs none.
- The [AMReX](https://github.com/AMReX-Codes/amrex) source, cloned next to COMPAS. It is compiled
  as part of each case, not built or installed separately.
- Python 3 with numpy, matplotlib, scipy and yt, only for the analysis and plotting scripts
  (see [Python environment](#python-environment)).
- For viewing results, ParaView or VisIt, or yt in Python. All are optional.

## Getting COMPAS and AMReX

Clone COMPAS and AMReX side by side in one workspace:

```bash
mkdir compas-workspace && cd compas-workspace
git clone https://github.com/scfpl-umich/COMPAS.git
git clone --branch 26.01 https://github.com/AMReX-Codes/amrex.git
```

```
compas-workspace/
├── amrex/
└── COMPAS/
    ├── source/        driver and model classes
    ├── include/       physics and operators, header-only
    ├── exec/
    │   ├── Make.exec  build rules shared by the cases
    │   ├── _Tests/    small cases run by the test suite
    │   └── ...        larger application cases
    └── scripts/       test runner and Python tools
```

Each case finds COMPAS from its own location and looks for AMReX in `../amrex` relative to the
COMPAS directory; if AMReX is elsewhere, add `AMREX_HOME=/path/to/amrex` to the `make` line.
Building a case kept outside the repository is described in
[Setting up a case](../user-guide/cases.md#starting-a-new-case).

### AMReX versions

COMPAS is tested with AMReX 26.01, the release the clone above checks out. The AMReX version of a
build is printed at startup and recorded in the configuration log of every run
([What a run writes](quickstart.md#what-a-run-writes)).

## Building

Build from a case directory with GNU Make:

```bash
cd COMPAS/exec/_Tests/Advection-5Eq
make -j4
```

The build ends with `SUCCESS`. The first build of a case compiles COMPAS and AMReX and takes a
minute or two on a laptop; later builds recompile only what changed. Object files go in
`tmp_build_dir/` inside the case, and `make realclean` removes them and every executable in the
case directory.

The model and the backends are set in the case's `GNUmakefile`. Variables such as `DIM`, `USE_MPI`
and `DEBUG` can also be set on the `make` line, which overrides the file:

```bash
make -j4 USE_MPI=FALSE     # serial build
make -j4 DEBUG=TRUE        # assertions and array bounds checks
```

Every option is listed in [Compile-time options](../user-guide/compile-options.md).

### The executable name

The executable is written to the case directory, and AMReX names it after the build:

```
main2d.gnu.MPI.ex
    │   │   └── parallel backends
    │   └────── compiler family (COMP)
    └────────── dimension (DIM)
```

A few examples:

| Build | Executable |
|---|---|
| 2D, MPI (the default of most cases) | `main2d.gnu.MPI.ex` |
| 2D, serial (`USE_MPI=FALSE`) | `main2d.gnu.ex` |
| 2D, MPI, debug (`DEBUG=TRUE`) | `main2d.gnu.DEBUG.MPI.ex` |
| 3D, MPI, OpenMP and CUDA | `main3d.gnu.MPI.OMP.CUDA.ex` |
| 3D, MPI and HIP (`USE_HIP=TRUE`) | `main3d.hip.MPI.HIP.ex` |

Builds with different options therefore sit side by side in the case directory. A serial build
runs without `mpirun`.

## Python environment

The analysis in these pages uses [yt](https://yt-project.org/) to read AMReX plotfiles, with
numpy, scipy, matplotlib and contourpy, kept in one virtual environment,
`scripts/compas_python_env`, which git ignores. Run outside a virtual environment,
`scripts/plot.py`, `scripts/convergence.py`, `scripts/postprocess.py` and the case analysis
scripts that import it create the environment and install the packages on first use. To create it
yourself:

```bash
cd COMPAS
python3 -m venv scripts/compas_python_env
scripts/compas_python_env/bin/pip install -r scripts/requirements.txt
source scripts/compas_python_env/bin/activate      # in each new shell
```

Check that it works with

```bash
python -c "import yt, scipy; print(yt.__version__, scipy.__version__)"
```

## Platform notes

### macOS

Install the command-line tools and MPI with Homebrew:

```bash
xcode-select --install
brew install mpich python coreutils
```

`coreutils` provides `gtimeout`, which the test runner uses to limit the time of each run.

The cases set `COMP = gnu`; on macOS the MPI wrapper `mpicxx` calls Apple Clang, and the build
works without changes.

### Linux

On Ubuntu or Debian:

```bash
sudo apt install build-essential mpich python3-venv
```

### Clusters

Load a compiler, MPI and Python through the site's module system. The module names vary by site:

```bash
module load gcc openmpi python
```

AMReX distributes the patches of the mesh across MPI ranks, so choose `amr.max_grid_size` to give
each rank several patches. For OpenMP, build with `USE_OMP=TRUE` and set `OMP_NUM_THREADS`.
`exec/BubbleFreeSurface-WaterHeAir/run.slurm` is a Slurm job script to start from: it runs the
first `main*.ex` in the case directory on `SLURM_NTASKS` ranks with `srun`.

### GPUs

Set one of `USE_CUDA`, `USE_HIP` or `USE_SYCL` to `TRUE` and rebuild. The kernels are the same
on CPUs and GPUs, and one MPI rank per GPU is typical. For example, on an AMD machine:

```bash
make -j16 USE_HIP=TRUE USE_MPI=TRUE
```

A HIP build uses the `hip` compiler family and a SYCL build the `sycl` family, whatever `COMP`
is, and the executable name shows it.
