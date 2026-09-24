KelvinHelmholtz-5Eq: mode-amplitude histories and linear growth rates
=====================================================================

This directory holds the in-situ history of the seeded mode for each run of the case and the
growth rates fitted from them, which are the data of the growth-rate figure (kh_growth.png).
Every run uses the 2D MPI build of exec/KelvinHelmholtz-5Eq (make -j3, AMReX 26.01) on 3 MPI
ranks unless its entry says otherwise, with prob/inputs and the command-line overrides given
below. Commands are run from exec/KelvinHelmholtz-5Eq, with this directory as data/. Each run
writes plot/<case_name>/mode.csv, and the mode_*.csv files here are those files, renamed after
the run.

Columns of mode_*.csv, one line every run.user_output_int = 4 coarse steps, starting at t = 0
(see prob/ProblemICBC.H; nondimensional units, U = delta = rho = 1):
  t        time, in delta/U
  step     coarse step
  S, C     int v sin(kx x) dA and int v cos(kx x) dA over the domain, kx = 2 pi/Lx
  A        amplitude of the seeded Fourier mode of v, integrated over y, 2 sqrt(S^2 + C^2)/Lx
  Ey       y kinetic energy, int rho v^2/2 dA
  Ey_max   largest rho v^2/2
  theta    momentum thickness of the layer, int (U^2 - u^2)/(4 U^2) dA / Lx


Mach sweep, finest cells about delta/16 (amr.max_level = 2)
-----------------------------------------------------------
mode_M0.4.csv   M = 0.4, delta/dx = 16.14, t = 0 to 100 (910 lines). prob/inputs is set for this run:
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs

mode_M0.2.csv   M = 0.2, delta/dx = 16.27, t = 0 to 55 (841 lines):
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.2 prob.Mach=0.2 \
        stop_time=55 amr.t_write=1000 \
        geometry.prob_lo="0.0 -19.665682 -1.0" geometry.prob_hi="14.749261 19.665682 1.0" amr.n_cell="60 160 8"

mode_M0.6.csv   M = 0.6, delta/dx = 16.02, t = 0 to 77.9 (522 lines):
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.6 prob.Mach=0.6 \
        stop_time=78 amr.t_write=1000 \
        geometry.prob_lo="0.0 -19.978332 -1.0" geometry.prob_hi="16.981582 19.978332 1.0" amr.n_cell="68 160 8"

mode_M0.8.csv   M = 0.8, delta/dx = 15.99, t = 0 to 120 (677 lines):
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.8 prob.Mach=0.8 \
        stop_time=120 amr.t_write=1000 \
        geometry.prob_lo="0.0 -40.036227 -1.0" geometry.prob_hi="22.520377 40.036227 1.0" amr.n_cell="90 320 8"


Resolution check at M = 0.4, finest cells about delta/32 (amr.max_level = 3)
----------------------------------------------------------------------------
mode_M0.4-d32.csv   M = 0.4, delta/dx = 32.29, t = 0 to 55 (487 lines):
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.4-d32 amr.max_level=3 \
        prob.y_ref="8 6 5" stop_time=55 amr.t_write=1000


Sharp tracer at M = 0.4, without and with THINC
-----------------------------------------------
Both set the tanh width of the tracer to one finest cell, prob.tracer_width = dx/delta
= 15.36231/248 = 0.0619448, and write plotfiles at the four snapshot times.

mode_M0.4-sharp.csv   WENO5, delta/dx = 16.14, t = 0 to 95 (862 lines):
    mpirun -n 3 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.4-sharp \
        prob.tracer_width=0.0619448 stop_time=95 amr.t_write="65 75 85 95"

mode_M0.4-sharp-THINC.csv   WENO5 and THINC, delta/dx = 16.14, t = 0 to 95 (862 lines), 1 MPI rank:
    mpirun -n 1 ./main2d.gnu.MPI.ex prob/inputs amr.case_name=./plot/KH-M0.4-sharp-THINC \
        prob.tracer_width=0.0619448 stop_time=95 amr.t_write="65 75 85 95" FiniteVolume.THINC=1 \
        FiniteVolume.THINC.Multi=1 FiniteVolume.THINC.FaceState=1 FiniteVolume.THINC.PhaseDensity=1


Growth rates
------------
Columns: label; M; delta_over_dx, delta over the finest cell size (4 n_x/Lx with two levels,
8 n_x/Lx with three); rate, the growth rate alpha c_i in U/delta; rate_err, its uncertainty; t0 and
t1, the fit window (see analysis/growth.py).

growth_rates.csv   the Mach sweep and the resolution check, and the growth-rate figure:
    python3 analysis/growth.py 0.2:0.2:16.27:data/mode_M0.2.csv 0.4:0.4:16.14:data/mode_M0.4.csv \
        0.6:0.6:16.02:data/mode_M0.6.csv 0.8:0.8:15.99:data/mode_M0.8.csv \
        0.4-d32:0.4:32.29:data/mode_M0.4-d32.csv --out data/growth_rates.csv --figure kh_growth.png

growth_rates_sharp.csv   the sharp-tracer pair:
    python3 analysis/growth.py 0.4-sharp:0.4:16.14:data/mode_M0.4-sharp.csv \
        0.4-sharp-THINC:0.4:16.14:data/mode_M0.4-sharp-THINC.csv --out data/growth_rates_sharp.csv
