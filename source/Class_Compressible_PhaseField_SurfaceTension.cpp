// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Surface tension (SURFACE_TENSION = true): capillary time step, capillary energy of the reflux and
// initial pressure projection. The face fluxes and the cell terms are in compute_dUdt_FV
// (Class_Compressible_PhaseField.cpp) and include/Physics_SurfaceTension.

#include <COMPAS.H>
#include <Kernels.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Gpu.H>
#include <LinearAlgebra.H>
#include <Tools.H>
#include <Operators_FV.H>
#include <Operators_FD.H>
#include <Class_Compressible_PhaseField.H>

#include <Macros.H>
#if (PHYSICS==FIVEEQS || PHYSICS==FIVEEQS_NPHASE || PHYSICS==SIXEQS || PHYSICS==SIXEQS_IE_NPHASE)
#include <IncludePhysics.H>
#include <Parm.H>
#include <Operators_LinearSystem.H>
#include <Physics_PhaseField.H>
#include <Operators_PhaseField.H>
#endif

#if (SURFACE_TENSION == true)
#include <Physics_SurfaceTension.H>
#endif

#include <limits>

using namespace amrex;


void
Compressible_PhaseField::PostInitFromScratch ()
{
#if (SURFACE_TENSION == true)
    if (h_parm->Physics_Parm.SurfaceTension_Parm.init_pressure == 1){
        SurfaceTension_InitialPressure();
    }
#endif
}


// Capillary time step (Brackbill, Kothe and Zemach 1992), dt <= sqrt(rho_bar h^3/(2 pi sigma_max)),
// with rho_bar the smallest mixture density of the interface cells (1e-3 < alpha_k < 1 - 1e-3 for
// some phase) of this level and h its smallest cell size. The value is local to this rank;
// ComputeDt takes the minimum over the ranks. No interface cell or sigma_max = 0: no limit.
Real
Compressible_PhaseField::SurfaceTension_EstTimeStep (int lev)
{
    Real dt_cap = std::numeric_limits<Real>::max();
#if (SURFACE_TENSION == true)
    SurfaceTension_Parameter const& ST = h_parm->Physics_Parm.SurfaceTension_Parm;
    if (!(ST.sigma_max > 0.0)){
        return dt_cap;
    }

    const Real Alpha_Interface = 1.0e-3;
    const Real Huge = std::numeric_limits<Real>::max();

    MultiFab const& U = dof_new[lev];
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(U,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& UArray = U.const_array(mfi);
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            bool Interface = false;
            for (int i_Phase = 0; i_Phase < ST_NPHASE; i_Phase++){
                Real Alpha = SurfaceTension_VolumeFraction(i, j, k, i_Phase, UArray);
                if (Alpha > Alpha_Interface && Alpha < 1.0 - Alpha_Interface){
                    Interface = true;
                }
            }
            return { Interface ? SurfaceTension_Density(i, j, k, UArray) : Huge };
        });
    }
    Real rho_bar = amrex::get<0>(reduce_data.value());

    if (rho_bar < Huge){
        const Real* dx = geom[lev].CellSize();
        Real h = dx[0];
        for (int idim = 1; idim < AMREX_SPACEDIM; idim++){
            h = amrex::min(h, dx[idim]);
        }
        dt_cap = std::sqrt(rho_bar*h*h*h/(2.0*PI*ST.sigma_max));
    }
#else
    amrex::ignore_unused(lev);
#endif
    return dt_cap;
}


// The reflux corrects the momentum of the coarse cells next to a finer level by dt Delta[div T]_h,
// the mismatch of the coarse and fine capillary fluxes. With energy_form = work the energy gets the
// matching work u.(dt Delta[div T]_h), with the velocity of the state before the reflux, as
// FVM_SurfaceIntegral_NC_Capillary does in the cell update. With energy_form = conservative the
// energy flux is refluxed with the conservative fluxes; SIXEQS then moves the share of phase 2,
// Y_2 dt Delta[div(u.T)]_h, from energy 1 to energy 2.
void
Compressible_PhaseField::SurfaceTension_RefluxSourceTerms (MultiFab & new_dof,
                                                           MultiFab const & state_before,
                                                           int lev,
                                                           Vector<std::unique_ptr<FluxRegister>> & fr)
{
#if (SURFACE_TENSION == true) && (NONCONSERVATIVE == true)
#if (PHYSICS != SIXEQS)
    if (h_parm->Physics_Parm.SurfaceTension_Parm.energy_form != CapillaryEnergy::Work){
        return;
    }
#endif
    constexpr int ncomp = NC_TERMS - INDEX_NC_CapillaryStress;

    MultiFab xi(grids[lev], dmap[lev], ncomp, 0);
    xi.setVal(Real(0.0));
    fr[lev+1]->Reflux(xi, 1.0, INDEX_NC_CapillaryStress, 0, ncomp, geom[lev]);

    Parm const* lparm = d_parm;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(new_dof,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& StateArray = state_before.const_array(mfi);
        Array4<Real const> const& xiArray    = xi.const_array(mfi);
        Array4<Real      > const& dofArray   = new_dof.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // xi = -dt Delta[div F]_h for each NC component
            Array<Real,AMREX_SPACEDIM> Div_Stress;
            for (int b = 0; b < AMREX_SPACEDIM; b++){
                Div_Stress[b] = -xiArray(i,j,k,b);
            }
            Real Div_Work = 0.0;
#if (PHYSICS == SIXEQS)
            Div_Work = -xiArray(i,j,k,AMREX_SPACEDIM);
#endif
            SurfaceTension_Energy_D(i, j, k,
                                    dofArray, StateArray,
                                    Div_Stress, Div_Work,
                                    lparm->Physics_Parm.SurfaceTension_Parm);
        });
    }
#else
    amrex::ignore_unused(new_dof, state_before, lev, fr);
#endif
}


// Volume fractions read by the capillary stress, on the valid cells grown by ngrow: the stored
// alpha_k of U (ST_NALPHA components from INDEX_VolumeFraction1), then SurfaceTension.smoothing
// passes of the 1-2-1 filter, each on one ghost layer less. Without smoothing Alpha is a copy, and
// the stress is the one of the stored volume fractions.
void
Compressible_PhaseField::SurfaceTension_VolumeFractions (MultiFab const & U,
                                                         MultiFab & Alpha,
                                                         int ngrow)
{
#if (SURFACE_TENSION == true)
    const int n_Smooth = h_parm->Physics_Parm.SurfaceTension_Parm.smoothing;
    const int ng = ngrow + n_Smooth;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(U.nGrow() >= ng && Alpha.nGrow() >= ngrow,
                                     "SurfaceTension_VolumeFractions: not enough ghost cells");

    if (n_Smooth == 0){
        MultiFab::Copy(Alpha, U, INDEX_VolumeFraction1, 0, ST_NALPHA, ngrow);
        return;
    }

    MultiFab A0(U.boxArray(), U.DistributionMap(), ST_NALPHA, ng);
    MultiFab A1(U.boxArray(), U.DistributionMap(), ST_NALPHA, ng);
    MultiFab::Copy(A0, U, INDEX_VolumeFraction1, 0, ST_NALPHA, ng);
    for (int pass = 1; pass <= n_Smooth; pass++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(A1,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.growntilebox(ng - pass);
            Array4<Real const> const& In  = A0.const_array(mfi);
            Array4<Real      > const& Out = A1.array(mfi);
            amrex::ParallelFor(bx, ST_NALPHA,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Out(i,j,k,n) = SurfaceTension_Smooth_K(i, j, k, n, In);
            });
        }
        std::swap(A0, A1);
    }
    MultiFab::Copy(Alpha, A0, 0, 0, ST_NALPHA, ngrow);
#else
    amrex::ignore_unused(U, Alpha, ngrow);
#endif
}


// =======================================
// Initial pressure projection (SurfaceTension.init_pressure = 1, fresh starts only), on all levels:
//  1. a ghosted copy of the state (two ghost cells);
//  2. the mixture pressure p0 and the cell capillary force f = [div T]_h, on the cells and the
//     first ghost layer;
//  3. face residuals R_f = beta_f ((f_L + f_R)/2 . e_d - (p0_R - p0_L)/dx_d), with
//     beta_f = 2/(rho_L + rho_R) (init_pressure_weight = density) or 1 (unit), and R_f = 0 on the
//     non-periodic domain faces;
//  4. R and beta averaged down to the coarse faces;
//  5. div(beta grad dp) = div(R) with the composite MLMG solve of the Phase-Field step (Neumann,
//     periodic where the domain is periodic), converged or aborting;
//  6. the volume-weighted mean of dp removed;
//  7. p = p0 + dp;
//  8. abort if p + pinf_k <= 0 for a present phase;
//  9. the energy rebuilt from p with alpha_k, alpha_k rho_k and u fixed (SurfaceTension_SetPressure_K);
// 10. average down, and dof_old = dof_new.
// The projection removes the gradient part of the mismatch between the discrete capillary force
// and the discrete pressure gradient; an exact discrete rest is not possible with the Riemann fluxes.
// =======================================
void
Compressible_PhaseField::SurfaceTension_InitialPressure ()
{
#if (SURFACE_TENSION == true)
    BL_PROFILE("Compressible_PhaseField::SurfaceTension_InitialPressure()");

    SurfaceTension_Parameter const& ST = h_parm->Physics_Parm.SurfaceTension_Parm;
    const bool Density_Weight = (ST.init_pressure_weight == CapillaryWeight::Density);
    const int nlevs = finest_level + 1;
    Parm const* lparm = d_parm;

    Vector<MultiFab> P0(nlevs);
    Vector<MultiFab> Rhs(nlevs);
    Vector<MultiFab> dP(nlevs);
    Vector<MultiFab> acoef(nlevs);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> Residual(nlevs);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef(nlevs);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> Flux_dP(nlevs);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> Grad_dP(nlevs);

    for (int lev = 0; lev < nlevs; lev++)
    {
        P0[lev].define(grids[lev], dmap[lev], 1, 1);
        Rhs[lev].define(grids[lev], dmap[lev], 1, 0);
        dP[lev].define(grids[lev], dmap[lev], 1, 1);
        dP[lev].setVal(0.0);
        acoef[lev].define(grids[lev], dmap[lev], 1, 0);
        acoef[lev].setVal(0.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
            BoxArray ba = grids[lev];
            ba.surroundingNodes(idim);
            Residual[lev][idim].define(ba, dmap[lev], 1, 0);
            bcoef[lev][idim].define(ba, dmap[lev], 1, 0);
            Flux_dP[lev][idim].define(ba, dmap[lev], 1, 0);
            Grad_dP[lev][idim].define(ba, dmap[lev], 1, 0);
        }

        // 1. Ghosted state: the stress on the faces of the first ghost layer reads the second one
        //    (and the smoothing one more layer per pass)
        MultiFab Sborder(grids[lev], dmap[lev], NSTATE, 2 + ST.smoothing);
        FillPatch(Sborder, lev, t_new[lev],
                  dof_new, t_new,
                  dof_old, t_old,
                  0, 0, NSTATE);
        MultiFab Alpha_ST(grids[lev], dmap[lev], ST_NALPHA, 2);
        SurfaceTension_VolumeFractions(Sborder, Alpha_ST, 2);

        MultiFab Rho(grids[lev], dmap[lev], 1, 1);
        MultiFab Force(grids[lev], dmap[lev], AMREX_SPACEDIM, 1);

        const Real dx = geom[lev].CellSize(0);
#if (AMREX_SPACEDIM > 1)
        const Real dy = geom[lev].CellSize(1);
#else
        const Real dy = 1.0;
#endif
#if (AMREX_SPACEDIM > 2)
        const Real dz = geom[lev].CellSize(2);
#else
        const Real dz = 1.0;
#endif
        amrex::Array<amrex::Real, AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx,dy,dz)};
        Box const& domain = geom[lev].Domain();

        for (MFIter mfi(P0[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.validbox();
            const Box  gbx = amrex::grow(bx, 1);

            Array4<Real const> const& U        = Sborder.const_array(mfi);
            Array4<Real const> const& A        = Alpha_ST.const_array(mfi);
            Array4<Real      > const& P0Array  = P0[lev].array(mfi);
            Array4<Real      > const& RhoArray = Rho.array(mfi);
            Array4<Real      > const& FArray   = Force.array(mfi);

            // 2. p0, rho and f = [div T]_h on the cells and the first ghost layer
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                P0Array(i,j,k)  = SurfaceTension_Pressure_K(i, j, k, U, *lparm);
                RhoArray(i,j,k) = SurfaceTension_Density(i, j, k, U);
                for (int b = 0; b < AMREX_SPACEDIM; b++){
                    FArray(i,j,k,b) = 0.0;
                }
                for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
                    Array<Real,AMREX_SPACEDIM> Stress_m;
                    Array<Real,AMREX_SPACEDIM> Stress_p;
                    FDM_CapillaryStress_K(i, j, k, Stress_m, A, AMREX_D_DECL(dx, dy, dz), idim,
                                          lparm->Physics_Parm.SurfaceTension_Parm);
                    FDM_CapillaryStress_K(i + (idim == 0), j + (idim == 1), k + (idim == 2),
                                          Stress_p, A, AMREX_D_DECL(dx, dy, dz), idim,
                                          lparm->Physics_Parm.SurfaceTension_Parm);
                    for (int b = 0; b < AMREX_SPACEDIM; b++){
                        FArray(i,j,k,b) += (Stress_p[b] - Stress_m[b])/dX[idim];
                    }
                }
            });

            // 3. Face residuals and face weights
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
            {
                const Box fbx = surroundingNodes(bx, idim);
                Array4<Real> const& RArray = Residual[lev][idim].array(mfi);
                Array4<Real> const& BArray = bcoef[lev][idim].array(mfi);
                const bool Periodic = geom[lev].isPeriodic(idim);
                const int  Face_lo  = domain.smallEnd(idim);
                const int  Face_hi  = domain.bigEnd(idim) + 1;
                const Real h        = dX[idim];
                amrex::ParallelFor(fbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    int im = (idim == 0) ? i-1 : i;
                    int jm = (idim == 1) ? j-1 : j;
                    int km = (idim == 2) ? k-1 : k;
                    Real Beta = Density_Weight ? 2.0/(RhoArray(im,jm,km) + RhoArray(i,j,k)) : 1.0;
                    BArray(i,j,k) = Beta;
                    int Face = (idim == 0) ? i : ((idim == 1) ? j : k);
                    bool Boundary = !Periodic && (Face == Face_lo || Face == Face_hi);
                    RArray(i,j,k) = Boundary ? 0.0
                                  : Beta*( 0.5*(FArray(im,jm,km,idim) + FArray(i,j,k,idim))
                                         - (P0Array(i,j,k) - P0Array(im,jm,km))/h );
                });
            }
        }
    }

    // 4. Average the face residuals and weights down to the coarse faces
    AverageDownFaces(Residual);
    AverageDownFaces(bcoef);

    // Right-hand side div(R)
    for (int lev = 0; lev < nlevs; lev++)
    {
        const Real* dx = geom[lev].CellSize();
        amrex::Array<amrex::Real, AMREX_SPACEDIM> dX;
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++){
            dX[idim] = dx[idim];
        }
        for (MFIter mfi(Rhs[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<Real> const& RhsArray = Rhs[lev].array(mfi);
            AMREX_D_TERM(Array4<Real const> const& Rx = Residual[lev][0].const_array(mfi);,
                         Array4<Real const> const& Ry = Residual[lev][1].const_array(mfi);,
                         Array4<Real const> const& Rz = Residual[lev][2].const_array(mfi));
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                RhsArray(i,j,k) = AMREX_D_TERM( (Rx(i+1,j,k) - Rx(i,j,k))/dX[0],
                                              + (Ry(i,j+1,k) - Ry(i,j,k))/dX[1],
                                              + (Rz(i,j,k+1) - Rz(i,j,k))/dX[2] );
            });
        }
    }

    // 5. Composite solve of div(beta grad dp) = div(R): a = 0, b = -1
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> LinOpBC_lo;
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> LinOpBC_hi;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        if ( geom[0].isPeriodic(idim) ){
            LinOpBC_lo[idim] = LinOpBCType::Periodic;
            LinOpBC_hi[idim] = LinOpBCType::Periodic;
        }
        else{
            LinOpBC_lo[idim] = LinOpBCType::Neumann;
            LinOpBC_hi[idim] = LinOpBCType::Neumann;
        }
    }
    LinearSystem_Parameter LinearSystem_Parm = h_parm->LinearSystem_Parm;
    LinearSystem_Parm.composite_solve = true;
    LinearSystem_Parm.fixed_iter      = 0;     // MLMG aborts if it does not converge
    LinearSystem_Parm.max_iter        = 1000;
    LinearSystem_Parm.max_fmg_iter    = 0;
    amrex::Real ascalar =  0.0;
    amrex::Real bscalar = -1.0;

    Print() << "SurfaceTension.init_pressure: projecting the initial pressure on " << nlevs << " level(s)\n";
    SolveABecLaplacianAllLevels_HomBC(dP, Flux_dP, Grad_dP,
                                      Geom(0, finest_level), grids, dmap, ref_ratio,
                                      ascalar, bscalar, acoef, bcoef, Rhs,
                                      LinOpBC_lo, LinOpBC_hi,
                                      0, nullptr, 2,
                                      LinearSystem_Parm);

    // 6. Remove the volume-weighted mean (level 0 holds the composite average after averaging down)
    for (int lev = finest_level-1; lev >= 0; --lev){
        amrex::average_down(dP[lev+1], dP[lev], geom[lev+1], geom[lev], 0, 1, refRatio(lev));
    }
    const Real Mean = dP[0].sum(0)/Real(geom[0].Domain().numPts());
    for (int lev = 0; lev < nlevs; lev++){
        dP[lev].plus(-Mean, 0, 1, 0);
    }
    const Real dP_max = dP[0].norm0(0);

    // 7-9. p = p0 + dp and the energy rebuilt from it
    int n_Fail = 0;
    for (int lev = 0; lev < nlevs; lev++)
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(dof_new[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<Real      > const& UArray  = dof_new[lev].array(mfi);
            Array4<Real const> const& P0Array = P0[lev].const_array(mfi);
            Array4<Real const> const& dPArray = dP[lev].const_array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real p = P0Array(i,j,k) + dPArray(i,j,k);
                return { SurfaceTension_SetPressure_K(i, j, k, UArray, p, *lparm) };
            });
        }
        n_Fail += amrex::get<0>(reduce_data.value());
    }
    ParallelDescriptor::ReduceIntSum(n_Fail);
    if (n_Fail > 0){
        amrex::Abort("SurfaceTension.init_pressure = 1: the projected pressure p gives p + pinf_k <= 0 for a present "
                     "phase in " + std::to_string(n_Fail) + " cells, where no state exists. Raise the ambient "
                     "pressure, or start without the projection.");
    }

    // 10. Consistent levels, and the old state equal to the new one
    AverageDown(dof_new);
    for (int lev = 0; lev < nlevs; lev++){
        MultiFab::Copy(dof_old[lev], dof_new[lev], 0, 0, NSTATE, 0);
    }

    Print() << "SurfaceTension.init_pressure: max |dp| = " << dP_max
            << " (mean " << Mean << " removed)\n";
#endif
}
