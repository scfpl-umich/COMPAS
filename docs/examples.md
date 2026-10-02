# Examples

The application cases are in `exec/`, grouped below by model. Each is set up and tested in 2D.
The header of its `prob/inputs` describes the problem and its units, and the inputs set the methods,
so the commands on these pages show only what a run changes. A case builds and runs like the
[first simulation](getting-started/quickstart.md), from its own directory:

```bash
cd exec/DropletSplashing-5Eq
make -j4
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs
```

Every case except `IsentropicVortex-5Eq` writes in-situ diagnostics to a CSV file in its plot
directory through its `UserOutputFunction`
([Setting up a case](user-guide/cases.md#in-situ-diagnostics)).

## Cases compared with a reference

These cases are compared with exact solutions, theory or experiments on the
[Verification and validation](verification.md) page:

| Case | Model | Problem | Also shows |
|---|---|---|---|
| [`DoubleMachReflection-5Eq`](verification.md#double-mach-reflection) | `FIVEEQS` | Mach 10 shock reflecting off a wedge | a `BoundaryFill_User` face that follows the exact moving shock |
| [`GuderleyImplosion-5Eq`](verification.md#guderley-converging-shock) | `FIVEEQS` | Converging shock on a perturbed heavy gas | the similarity solution as initial condition and exact inflow |
| [`VortexWindup-5Eq`](verification.md#vortex-windup) | `FIVEEQS` | Interface wound up by a steady vortex | THINC in a sheared flow, with an exact solution |
| [`ThreeMaterialRM-5Eq-N`](verification.md#shock-through-two-interfaces) | `FIVEEQS_NPHASE`, 3 phases | Shock through two successive interfaces | SI units; `make NPHASE=4` adds an auxiliary material for a reduction-consistency test |
| [`ShockBubble-AirHelium-5Eq`](verification.md#shock-and-helium-cylinder) | `FIVEEQS` | Shock and helium cylinder | six levels of refinement on the interface and the shocks |
| [`ShockRefraction-AirWater-5Eq`](verification.md#shock-refraction-at-an-air-water-interface) | `FIVEEQS` | Shock refraction at a water surface | a stiffened gas, a sharp initial interface and a post-shock inflow |
| [`KelvinHelmholtz-5Eq`](verification.md#compressible-kelvin-helmholtz-instability) | `FIVEEQS` | Compressible shear layer | refined bands from `spatio_temporal_tag`; the run histories in `data/` |
| [`WaterHammer-6Eq`](verification.md#planar-water-hammer) | `SIXEQS` | Water column stopped by a wall | the exact solution in `Parm::DynamicInit` |

The test cases in `exec/_Tests/` that appear there are listed under [Test cases](#test-cases).

## Two-phase five-equation model

### DropletSplashing-5Eq

`PHYSICS=FIVEEQS`, `DIFFUSION=true`, 2D. A heavy droplet of radius 0.25 falls under gravity
through a light gas into a pool of the same heavy material, in a closed box, with density ratio
100. The case shows the viscous five-equation model with gravity added through
`user_source_term`, and writes the area, centroid height and mean vertical velocity of the heavy
fluid above $y = R/2$ to `droplet.csv`. The figures come from two runs with one more level than the
inputs, 512 by 768 finest cells, with MUSCL-MC alone and with THINC:

```bash
cd exec/DropletSplashing-5Eq && make -j4
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs amr.max_level=3 amr.t_write_interval=0.125
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs amr.max_level=3 amr.t_write_interval=0.125 \
    FiniteVolume.THINC=1 amr.case_name=./plot/DropletSplashing-5Eq-THINC
```

![DropletSplashing-5Eq with MUSCL-MC and with THINC at three times](media/figures/dropletsplashing_thinc_snapshots.png)

*`exec/DropletSplashing-5Eq` at impact ($t = 2$), in the crater ($t = 2.5$) and late in the
splash ($t = 3$), with MUSCL-MC (top row) and MUSCL-MC + THINC (bottom row). Right half of each
panel: light-gas volume fraction $\alpha_1$. Left half: pressure schlieren $|\nabla P|$, the
interface $\alpha_1 = 0.5$ (black) and the outlines of AMR levels 1 to 3.*

By impact the droplet has flattened and the gas flowing around it has rolled its edges into thin
rims. It traps a thin gas film against the pool, merges with the pool and opens a crater, whose rim
throws up a splash sheet as the pool climbs the side walls. With THINC the droplet and pool
surfaces stay about two cells thick until impact and about three afterwards, against 6 to 12 cells
with MUSCL-MC alone. The difference is largest in the gas film, which MUSCL-MC spreads into a
diffuse band and THINC keeps sharp until it breaks into a chain of small gas pockets. With THINC
the splash is a tall sheet whose rim rolls up and sheds small droplets; without it, it forms
thicker, rolled-up fingers. THINC leaves small roll-ups on the sides of the falling droplet, where
the gas flows past, and the thinnest late sheets, below a cell or two, as faint smeared trails.

![Heavy-fluid area, centroid height and mean vertical velocity above y = R/2 in the two runs](media/figures/dropletsplashing_thinc_history.png)

*Heavy fluid above the cut $y = R/2$, with MUSCL-MC and with MUSCL-MC + THINC. (a) Area $V_2$,
normalized by the droplet area $\pi R^2$ (dashed); the rise after $t \approx 2.7$ is the pool
climbing the walls past the cut. (b) Centroid height $y_2$, with free fall from rest (dashed) and
the cut height (dash-dot). (c) Mass-weighted mean vertical velocity $v_2$, with free fall (dashed).
Dotted line: the free-fall contact time $t_c$.*

The diagnostics are nearly the same in both runs. Until impact the heavy fluid above the cut keeps
the area of the droplet to within 1 percent, and its centroid and velocity follow free fall,
trailing it slightly because of the drag of the gas. After impact its mean velocity turns upward
with the splash rim, peaking at 0.72 with MUSCL-MC and 0.59 with THINC.

### GasImpact-5Eq

`PHYSICS=FIVEEQS`, 2D. A dense, pressurized rectangular block moving at $u = 0.2$ strikes a
plate 15 times denser than the surrounding gas. The plate and the gas share one EOS and differ
only in density, and symmetry at $y = 1$ halves the block. The case shows THINC under a strong
density contrast. `GasImpact-5Eq-N` is the three-material version.

### CylindricalRM-5Eq

`PHYSICS=FIVEEQS`, 2D. A Mach 2 cylindrical shock in air converges on a heavy SF$_6$-like gas
whose interface is $R(\theta) = 0.5 + 0.02 \cos(8\theta)$, in one quarter of the plane with
symmetry on the axes. The case shows THINC in a converging flow, AMR on the interface and the
shocks, and in-situ diagnostics: it writes the equivalent radius of the heavy gas and the amplitude
of the perturbation mode to `interface.csv`, and the scripts in `analysis/` plot the contours, the
time series and the tracked amplitude. The COMPAS tutorial builds this case step by step, and the
version here is the finished one. `CylindricalRM-5Eq-N` is the three-material version.

## N-phase five-equation model

### CylindricalRM-5Eq-N

`PHYSICS=FIVEEQS_NPHASE`, `NPHASE=3`, 2D. The converging shock of `CylindricalRM-5Eq` meets two
interfaces: a Mach 2 cylindrical shock starts at $r = 0.92$ in air and converges on a heavy
SF$_6$-like shell with outer interface $r = 0.75 + 0.03 \cos(16\theta)$, around a light
helium-like core with interface $r = 0.5 + 0.02 \cos(12\theta)$. All three are ideal gases in
nondimensional units, and one quarter of the plane is simulated. The case shows the N-material
model with THINC on both interfaces and AMR that follows the interfaces and the shocks, and writes
the equivalent radii of the two interfaces and the half widths of their mixing zones to
`interfaces.csv`. The video on the [home page](index.md) and the figure use four levels above the
64 by 64 base grid, 1024 by 1024 cells at the finest level; the inputs use three:

```bash
cd exec/CylindricalRM-5Eq-N && make -j4
mpirun -n 8 ./main2d.gnu.MPI.ex prob/inputs amr.max_level=4 amr.t_write_interval=0.008
```

![The three-material implosion at five times](media/figures/cylrm3_sequence.png)

*`exec/CylindricalRM-5Eq-N` at $t = 0.2$, 0.3, 0.4, 0.56 and 1: air sand, SF$_6$ shell teal and
helium core dark blue, and in the lower half of each panel the pressure-gradient schlieren, the two
interfaces (black) and the outline of the region each AMR level covers. The panels zoom in as the
shell converges; each has its own scale bar, of length 0.1.*

By $t = 0.2$ to 0.3 the shock has crossed both interfaces; the shell perturbations grow into
spikes while the core stays nearly round. By $t = 0.4$ the shell has become long SF$_6$ jets with
rolled-up mushroom caps and small Kelvin-Helmholtz roll-ups along their stems, the helium core bulges
out between their roots, and the shocks focus inside the core. Near $t = 0.56$ the reflected shock
comes back out and re-shocks the debris, and the imploded region is smallest near $t = 0.57$. After
reshock the jets break up, and SF$_6$ and helium mix into a zone that grows back out to about
$r = 0.5$ by $t = 1$, with parcels down to a few cells across. THINC keeps the interfaces about two
cells thick through the jets and the roll-ups, and about two and a half in the mixing zone, and the
finest level follows the interfaces and the shocks, so the schlieren shows the transmitted and
reflected shocks inside the shell and the core.

### GasImpact-5Eq-N

`PHYSICS=FIVEEQS_NPHASE`, `NPHASE=3`, 2D. The block of `GasImpact-5Eq` strikes a plate of a
second material surrounded by a light gas, a third material, each with its own EOS. The case shows
the N-material model with THINC. The video compares the case as provided with the same inputs
without THINC, both with three levels of refinement, 1152 by 512 finest cells:

```bash
cd exec/GasImpact-5Eq-N && make -j4
mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/thinc amr.chk_int=-1
mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs FiniteVolume.THINC=0 amr.case_name=./plot/base amr.chk_int=-1
```

```{raw} html
<figure>
<video controls preload="metadata" playsinline muted loop width="1200" style="max-width:100%;height:auto"
       poster="_static/gasimpact_thinc_vs_base.png">
  <source src="_static/gasimpact_thinc_vs_base.mp4" type="video/mp4">
  <a href="_static/gasimpact_thinc_vs_base.mp4">Download the video (MP4)</a>
</video>
<figcaption>The block (dark blue) striking the plate (teal) in the gas (sand), from t = 0 to 10:
the THINC run above the symmetry plane y = 1 and the run with MUSCL-MC alone mirrored below it,
in colors blended by volume fraction.
<a href="_static/gasimpact_thinc_vs_base.mp4">Download the MP4</a>.</figcaption>
</figure>
```

At first the two halves are nearly the same: the corners of the block roll up, and weak circular
waves run through the gas. As the block pushes into the plate the runs differ behind it: without
THINC the block material rolls up into smooth vortices and mixes diffusely with the gas, and with
THINC the interfaces stay sharp and the block leaves a cloud of parcels a few cells across. The
outer shape of the block, its leading edge and the height of its main body agree closely. Late in
the run the block spreads into a thin cap wrapped in plate material, and the stem of the plate
bends back toward the axis, smooth without THINC and broken into sharp fragments with it.

### TriplePoint-5Eq-N

`PHYSICS=FIVEEQS_NPHASE`, `NPHASE=3`, 2D. Three gases at rest in a $7 \times 3$ box with
reflecting walls (Galera, Maire and Breil, J. Comput. Phys. 229, 2010): a high-pressure driver
($\rho = 1$, $P = 1$, $\gamma = 1.5$) fills $x < 1$, and to its right lie a dense gas below
$y = 1.5$ ($\rho = 1$, $P = 0.1$, $\gamma = 1.4$) and a light gas above ($\rho = 0.125$,
$P = 0.1$, $\gamma = 1.5$). The case shows the N-material model with sharp initial interfaces. The
inputs run it to $t = 5$ on 224 by 96 base cells with two levels, MUSCL-MC, HLLC and TVD-RK3, and
the second run adds THINC:

```bash
cd exec/TriplePoint-5Eq-N && make -j4
mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs
mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs FiniteVolume.THINC=1 amr.case_name=./plot/TriplePoint-THINC
```

![The three gases at t = 1, 3 and 5 without and with THINC](media/figures/triplepoint_sequence.png)

*`exec/TriplePoint-5Eq-N` at $t = 1$, 3 and 5, with MUSCL-MC (top row) and MUSCL-MC + THINC
(bottom row): dense gas dark blue, driver teal and light gas sand, blended by volume fraction, with
the outlines of AMR levels 1 and 2.*

The shock runs faster through the light gas than through the dense gas, so its front tilts, and
the shear along the dense-light interface puts vorticity at the triple point, where the dense gas
rolls up into a large spiral that wraps driver and light gas around it, as in the reference
results, and the shear layer behind it rolls up into a row of small Kelvin-Helmholtz billows. The
vortex and the fronts are in the same places with both schemes. With MUSCL-MC alone the interfaces
spread over nine to ten cells by $t = 5$, and the core of the roll-up is a diffuse blend of all
three gases. With THINC every interface stays under two cells thick, the billows stay sharp, and
the wound-up layers stay separate and break into thin filaments and small pockets of light gas and
driver.

## Six-equation models

### ShockBubble-WaterAir-6Eq

`PHYSICS=SIXEQS`, 2D. A Mach 1.2 shock in water (stiffened gas) hits a cylindrical air bubble
(ideal gas) of diameter 1, with the upper half simulated. The case shows the six-equation model
with instantaneous pressure relaxation, THINC and AMR at a liquid-gas interface.

### ShockDroplet-AirWater-6Eq

`PHYSICS=SIXEQS`, 2D. A Mach 4 shock in air hits a cylindrical water droplet of diameter 1, with
the upper half simulated. The case shows the six-equation model with pressure relaxation and AMR
on the droplet and the shocks. The figure uses THINC on a coarser grid than the inputs, 96 by 48
base cells with two levels, run to half the case's end time:

```bash
cd exec/ShockDroplet-AirWater-6Eq && make -j3
mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs FiniteVolume.THINC=1 stop_time=0.1 amr.n_cell="96 48 8" amr.max_level=2
```

![The shocked water column at four times, with THINC](media/figures/shockdroplet_snapshots.png)

*`exec/ShockDroplet-AirWater-6Eq` with THINC at $t u_0/D = 0.01$, 0.04, 0.07 and 0.10, mirrored
about $y = 0$. Upper half: air volume fraction $\alpha_1$. Lower half: pressure schlieren
$|\nabla P|$, the interface $\alpha_1 = 0.5$ and the patch boxes of AMR levels 1 and 2. The shock
moves in $+x$.*

A curved shock reflects from the upstream face of the column and becomes nearly planar by
$t = 0.1$, and the two branches of the incident shock meet in Mach stems downstream. The water
barely compresses, but the flow behind the shock flattens the column into a lens, with thin sheets
starting to peel off near its widest points by $t = 0.1$. THINC keeps the water-air interface one
to two cells thick the whole time, even on these coarse cells.

### ShockBubble-WaterAir-6Eq-Nphase

`PHYSICS=SIXEQS_IE_NPHASE`, `NPHASE=2`, 2D. The setup of `ShockBubble-WaterAir-6Eq` at Mach 1.1,
with the N-phase six-equation model in internal-energy form and instantaneous pressure relaxation.
The shock drives the upstream side of the bubble into a re-entrant jet that crosses it, and the air
is compressed to under one percent of its initial volume before it rebounds. The figures compare
the case as provided, with MUSCL, and the same inputs with THINC:

```bash
cd exec/ShockBubble-WaterAir-6Eq-Nphase && make -j4
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs
mpirun -n 4 ./main2d.gnu.MPI.ex prob/inputs FiniteVolume.THINC=1 \
    amr.case_name=./plot/ShockBubble-WaterAir-6Eq-Nphase-THINC
```

![ShockBubble-WaterAir-6Eq-Nphase with MUSCL and with THINC at four times through the collapse](media/figures/sbnphase_thinc_snapshots.png)

*`exec/ShockBubble-WaterAir-6Eq-Nphase` at $t u_0/D = 0.010$, 0.018, 0.020 and 0.024, with MUSCL
(top row) and MUSCL + THINC (bottom row). Upper half of each panel: air volume fraction
$\alpha_2$. Mirrored lower half: pressure schlieren $|\nabla P|$ with the interface
$\alpha_2 = 0.5$.*

![Air volume and axial interface positions with MUSCL and with THINC](media/figures/sbnphase_thinc_history.png)

*The two runs, with times in units of $D/u_0$, $u_0 = \sqrt{P_0/\rho_0}$. (a) Air volume $V$
over the simulated upper half, normalized by its initial value $V_0$; the inset shows it on a log
scale around the collapse, with the minima, 0.0063 (MUSCL) and 0.0059 (THINC), at
$t u_0/D = 0.020$ (open circles). (b) Positions on the axis of the upstream interface, which
becomes the jet, and of the downstream interface, until the jet reaches the downstream side at
$t_j$ (open circles). Dashed line: $t_a$, when the shock reaches the bubble.*

With THINC the bubble surface stays about two cells wide for the whole run, while with MUSCL alone
the upstream interface, which turns into the jet, spreads over many cells as the jet crosses the
bubble. After the collapse the remaining air is a single diffuse pocket trailing a spiral of partly
mixed fluid without THINC, and a few sharp pockets and small fragments with it. The air volume
follows nearly the same history in both runs.

### BubbleFreeSurface-WaterHeAir

`PHYSICS=SIXEQS_IE_NPHASE`, `NPHASE=3`, 2D. A gas bubble at 1/30 of the ambient pressure
collapses in water at 5 MPa, three bubble radii below a free surface with air above it. Half of
the bubble is simulated, with symmetry at $x = 0$. The case shows the N-phase six-equation model
with pressure relaxation on three materials and two interfaces of large density ratio. The same
inputs also run with `PHYSICS=FIVEEQS_NPHASE`, and `run.slurm` is a template for a cluster job.

## Test cases

The cases in `exec/_Tests/` are small. The [test suite](user-guide/testing.md) builds and runs
all of them, and each is also a starting point for a new case.

| Case | Model | Dim | Problem and what it tests |
|---|---|---|---|
| `Advection-5Eq` | `FIVEEQS` | 2 | Four shapes advected in a periodic box ([interface advection](verification.md#interface-advection)); base test of the model, with variants for every method |
| `Advection-5Eq-N` | `FIVEEQS_NPHASE`, 5 phases | 2 | Four shapes of four materials; base test of the model, with THINC |
| `Advection-6Eq` | `SIXEQS` | 2 | Four shapes; base test of the model, with pressure relaxation and ACDI |
| `Advection-6Eq-N` | `SIXEQS_IE_NPHASE`, 5 phases | 2 | Four shapes of four materials; base test of the model, with pressure relaxation and ACDI |
| `Advection-Viscous-5Eq-N` | `FIVEEQS_NPHASE`, 5 phases | 2 | Four shapes with viscosity and heat conduction in every material |
| `ShearDecay-6Eq` | `SIXEQS` | 2 | Decay of a shear wave in a mixture of two gases of different viscosities; viscous fluxes, the phase share of the dissipation, AMR reflux |
| `ShearDecay-6Eq-N` | `SIXEQS_IE_NPHASE`, 2 phases | 2 | As `ShearDecay-6Eq` |
| `HeatConduction-6Eq` | `SIXEQS` | 2 | Heat conduction across a contact between the two phases, against the erf profile; phase heat fluxes with pressure-temperature relaxation |
| `HeatConduction-6Eq-N` | `SIXEQS_IE_NPHASE`, 2 phases | 2 | The same contact, where no heat passes between the phases, and a uniform mixture against the erf profile |
| `Advection-Viscous-6Eq` | `SIXEQS` | 2 | A liquid drop advected in gas with viscosity and heat conduction; uniform pressure, velocity and phase temperatures |
| `Advection-Viscous-6Eq-N` | `SIXEQS_IE_NPHASE`, 3 phases | 2 | As `Advection-Viscous-6Eq`, with a third phase that stays exactly absent |
| `Advection-THINC-5Eq` | `FIVEEQS` | 2 | `Advection-5Eq` with THINC |
| `Advection-PF-5Eq` | `FIVEEQS` | 2 | `Advection-5Eq` with implicit CAC-Adv Phase-Field regularization |
| `COMPAS-STL-5Eq` | `FIVEEQS` | 2 | A solid read from an STL file falls into a pool; STL reader, gravity source term, viscosity, THINC, on CPUs |
| `RichtmyerMeshkov-5Eq` | `FIVEEQS` | 2 | Single-mode Richtmyer-Meshkov instability ([against experiment](verification.md#single-mode-richtmyer-meshkov-instability)); interface diagnostics with ghost cells, checkpoint and restart |
| `RichtmyerMeshkov-Multimode-5Eq` | `FIVEEQS` | 2 | Multimode Richtmyer-Meshkov instability; parameter arrays sized at run time and copied to the GPU |
| `RayleighTaylor-5Eq` | `FIVEEQS` | 2 | Single-mode Rayleigh-Taylor instability ([against theory](verification.md#single-mode-rayleigh-taylor-instability)); gravity source term, viscosity, THINC |
| `RayleighTaylor-3D-5Eq` | `FIVEEQS` | 3 | 3D Rayleigh-Taylor instability; 3D build, gravity source term, viscosity |
| `Jet-Inflow-5Eq` | `FIVEEQS` | 2 | Gas jet injected into a box with an outflow face; Robin user boundary function on an inflow face |
| `NonsphericalCollapse-6Eq` | `SIXEQS` | 3 | Gas bubble collapsing near a wall (Johnsen and Colonius 2009); pressure relaxation, THINC, user refinement tagging, no subcycling |
| `ShockVortex-5Eq` | `FIVEEQS` | 2 | Vortex through a stationary shock; quadrature WENO5 with AMR |
| `ViscousShockTube-5Eq` | `FIVEEQS` | 2 | Viscous shock tube of Daru and Tenaud; quadrature WENO5 with viscous fluxes and no-slip walls |
| `ViscousShockTube-6Eq` | `SIXEQS` | 2 | `ViscousShockTube-5Eq` with one gas in the six-equation model; WENO5 with viscous fluxes and no-slip walls |
| `ViscousShockTube-6Eq-N` | `SIXEQS_IE_NPHASE`, 2 phases | 2 | As `ViscousShockTube-6Eq` |
| `IsentropicVortex-5Eq` | `FIVEEQS` | 2 | Isentropic vortex over one period; quadrature WENO5, the [order-of-accuracy](verification.md#order-of-accuracy) case |
| `Sod-5Eq` | `FIVEEQS` | 3 | [Sod shock tube](verification.md#sod-shock-tube) with AMR and outflow boundaries; variants along $y$ and $z$, between walls, and the [water-air shock tube](verification.md#water-air-shock-tube) |
| `Sod-5Eq-N` | `FIVEEQS_NPHASE`, 2 phases | 3 | As `Sod-5Eq` |
| `Sod-6Eq` | `SIXEQS` | 3 | As `Sod-5Eq`, with pressure relaxation |
| `Sod-6Eq-N` | `SIXEQS_IE_NPHASE`, 2 phases | 3 | As `Sod-5Eq`, with pressure relaxation |
| `ShuOsher-5Eq` | `FIVEEQS` | 2 | [Shu-Osher](verification.md#shu-osher-problem) shock-entropy wave interaction with AMR; variants compare the reconstruction schemes |
| `NASG-5Eq` | `FIVEEQS` | 2 | Noble-Abel stiffened gas: air-water and water shock tubes with liquid water of Le Métayer and Saurel (2016), a water drop advected at uniform pressure and velocity, and $b = q = 0$ against the stiffened gas; conservation and covolume diagnostics in CSV |
| `NASG-5Eq-N` | `FIVEEQS_NPHASE`, 2 phases | 2 | As `NASG-5Eq` |
| `NASG-6Eq` | `SIXEQS` | 2 | As `NASG-5Eq`, with pressure and pressure-temperature relaxation |
| `NASG-6Eq-N` | `SIXEQS_IE_NPHASE`, 2 phases | 2 | As `NASG-5Eq`, with pressure relaxation |
