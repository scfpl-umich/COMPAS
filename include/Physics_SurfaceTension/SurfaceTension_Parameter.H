// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SurfaceTension_Parameter_H_
#define SurfaceTension_Parameter_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>
#if (SURFACE_TENSION == true)
#if (NONCONSERVATIVE != true)
#error "SURFACE_TENSION=true needs NONCONSERVATIVE=true: the capillary work and its reflux use the non-conservative face components."
#endif

#include<cmath>
#include<string>
#include<vector>

// =======================================
// Surface tension (compile with -DSURFACE_TENSION=true), read from ParmParse("SurfaceTension"):
//   sigma                 two-phase models: sigma_12; N-phase models: the N(N-1)/2 pairwise
//                         coefficients (1,2),(1,3),...,(1,N),(2,3),...,(N-1,N). Default 0.
//   energy_form           work (default): the energy gets u.[div T]_h, the work of the capillary
//                         force; conservative: the energy flux -(u_f.T.e_d)
//   init_pressure         1: project the initial pressure onto the discrete capillary force
//                         (fresh starts only); default 0
//   init_pressure_weight  density (default, beta = 1/rho) or unit (beta = 1), the weight of the
//                         projection
//   smoothing             passes of the 1-2-1 filter applied to the volume fractions before the
//                         stress is computed, 0 to 4; default 2 (a smoothed color function as in
//                         Brackbill, Kothe and Zemach 1992; Williams, Kothe and Puckett 1998)
//
// The capillary stress on a face with normal e_d is T.e_d = sum_k sigma_k (|g_k| e_d - g_kd g_k/|g_k|),
// with g_k the face gradient of alpha_k and the phase coefficients sigma_k defined by
// sigma_kl = sigma_k + sigma_l (FDM_Conservative2FluxCapillary_K.H).
// =======================================

// Number of phases that carry a capillary coefficient. The two-phase models define NPHASE = 1.
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
#define ST_NPHASE 2
#else
#define ST_NPHASE NPHASE
#endif

// Number of stored volume fractions the stress reads, from INDEX_VolumeFraction1 on: alpha_1 in the
// two-phase models, alpha_1 to alpha_{N-1} in the N-phase models (g_N = -sum of the others)
#define ST_NALPHA (ST_NPHASE - 1)

/** Non-conservative face components of surface tension, the last ones of NC_TERMS (Macros.H) */
enum NCIndex_SurfaceTension : int {
#if (PHYSICS == SIXEQS)
    INDEX_NC_CapillaryStress = NC_TERMS - AMREX_SPACEDIM - 1, /**< T.e_d, AMREX_SPACEDIM components */
    INDEX_NC_CapillaryWork   = NC_TERMS - 1                   /**< (u.T)_n, energy_form = conservative */
#else
    INDEX_NC_CapillaryStress = NC_TERMS - AMREX_SPACEDIM      /**< T.e_d, AMREX_SPACEDIM components */
#endif
};

enum class CapillaryEnergy : int {
    Work,
    Conservative
};

enum class CapillaryWeight : int {
    Density,
    Unit
};

struct SurfaceTension_Parameter
{
    amrex::Array<amrex::Real, ST_NPHASE*ST_NPHASE> sigma_pair{}; // sigma_kl at [k*ST_NPHASE + l], 0 for k = l
    amrex::Array<amrex::Real, ST_NPHASE> sigma{};                // phase coefficients, sigma_kl = sigma_k + sigma_l
    amrex::Real sigma_12  = 0.0;                                 // two-phase models: the stress uses sigma_12 directly
    amrex::Real sigma_max = 0.0;                                 // largest sigma_kl, for the capillary time step

    CapillaryEnergy energy_form   = CapillaryEnergy::Work;
    int init_pressure             = 0;
    CapillaryWeight init_pressure_weight = CapillaryWeight::Density;
    int smoothing                 = 2;

    // Pairwise coefficient sigma_kl (0-based phases)
    AMREX_GPU_HOST_DEVICE
    amrex::Real Pair (int k, int l) const { return sigma_pair[k*ST_NPHASE + l]; }

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("SurfaceTension");

        // Pairwise coefficients
        constexpr int n_pair = ST_NPHASE*(ST_NPHASE - 1)/2;
        std::vector<amrex::Real> sigma_list(n_pair, 0.0);
        if (pp.contains("sigma")){
            if (pp.countval("sigma") != n_pair){
                amrex::Abort("SurfaceTension.sigma has " + std::to_string(pp.countval("sigma")) + " values, but "
                             + std::to_string(ST_NPHASE) + " phases need " + std::to_string(n_pair)
                             + " pairwise values, ordered (1,2),(1,3),...,(1,N),(2,3),...,(N-1,N).");
            }
            pp.getarr("sigma", sigma_list);
        }
        int n = 0;
        for (int k = 0; k < ST_NPHASE; k++){
            sigma_pair[k*ST_NPHASE + k] = 0.0;
            for (int l = k+1; l < ST_NPHASE; l++){
                amrex::Real s = sigma_list[n];
                if (!(s >= 0.0)){
                    amrex::Abort("SurfaceTension.sigma must be >= 0, but value " + std::to_string(n+1)
                                 + " (phases " + std::to_string(k+1) + " and " + std::to_string(l+1) + ") is "
                                 + std::to_string(s) + ".");
                }
                sigma_pair[k*ST_NPHASE + l] = s;
                sigma_pair[l*ST_NPHASE + k] = s;
                sigma_max = amrex::max(sigma_max, s);
                n++;
            }
        }
        sigma_12 = sigma_pair[1];

        // Phase coefficients from sigma_kl = sigma_k + sigma_l: split equally for two phases,
        // unique for three, least squares for more (sigma_k = (r_k - S)/(N - 2) with
        // r_k = sum_{l != k} sigma_kl and S = sum_{k<l} sigma_kl/(N - 1); exact for N = 3)
        if (ST_NPHASE == 2){
            sigma[0] = 0.5*sigma_12;
            sigma[1] = 0.5*sigma_12;
        }
        else {
            amrex::Real S = 0.0;
            for (int k = 0; k < ST_NPHASE; k++){
                for (int l = k+1; l < ST_NPHASE; l++){
                    S += Pair(k,l);
                }
            }
            S /= amrex::Real(ST_NPHASE - 1);
            for (int k = 0; k < ST_NPHASE; k++){
                amrex::Real r = 0.0;
                for (int l = 0; l < ST_NPHASE; l++){
                    if (l != k) { r += Pair(k,l); }
                }
                sigma[k] = (r - S)/amrex::Real(ST_NPHASE - 2);
            }
            amrex::Real residual = 0.0;
            for (int k = 0; k < ST_NPHASE; k++){
                for (int l = k+1; l < ST_NPHASE; l++){
                    residual = amrex::max(residual, std::abs(sigma[k] + sigma[l] - Pair(k,l)));
                }
            }
            if (residual > 1.0e-12*sigma_max){
                amrex::Abort("SurfaceTension.sigma: the pairwise coefficients are not of the form sigma_kl = sigma_k + sigma_l "
                             "(least-squares residual " + std::to_string(residual) + "), which the capillary stress "
                             "sum_k sigma_k T(grad alpha_k) needs. Every set of three phases fits; for more phases "
                             "choose phase coefficients sigma_k and give their pairwise sums.");
            }
            for (int k = 0; k < ST_NPHASE; k++){
                if (sigma[k] < 0.0){
                    amrex::Print() << "Warning: SurfaceTension.sigma gives the negative phase coefficient sigma_"
                                   << k+1 << " = " << sigma[k] << " (the pairwise coefficients violate the triangle "
                                   << "inequality, so no Neumann triangle exists).\n";
                }
            }
        }

        std::string EnergyString = "work";
        pp.query("energy_form", EnergyString);
        if (EnergyString == "work"){
            energy_form = CapillaryEnergy::Work;
        }
        else if (EnergyString == "conservative"){
            energy_form = CapillaryEnergy::Conservative;
        }
        else {
            amrex::Abort("SurfaceTension.energy_form = '" + EnergyString + "' is not recognized. Valid values are "
                         "work, conservative.");
        }

        pp.query("init_pressure", init_pressure);
        if (init_pressure != 0 && init_pressure != 1){
            amrex::Abort("SurfaceTension.init_pressure = '" + std::to_string(init_pressure) + "' is not recognized. "
                         "Valid values are 0, 1.");
        }

        std::string WeightString = "density";
        pp.query("init_pressure_weight", WeightString);
        if (WeightString == "density"){
            init_pressure_weight = CapillaryWeight::Density;
        }
        else if (WeightString == "unit"){
            init_pressure_weight = CapillaryWeight::Unit;
        }
        else {
            amrex::Abort("SurfaceTension.init_pressure_weight = '" + WeightString + "' is not recognized. Valid values are "
                         "density, unit.");
        }

        // Smoothing of the volume fractions of the stress. The stress of a profile only one to three
        // cells wide (THINC, or a tanh with eta = h) is rough at the grid scale, and its spurious
        // currents grow until the run fails (THINC in every model, SIXEQS with relaxation on fine
        // grids), so the default is 2 passes.
        pp.query("smoothing", smoothing);
        if (smoothing < 0 || smoothing > 4){
            amrex::Abort("SurfaceTension.smoothing = " + std::to_string(smoothing) + " is not valid. Valid values "
                         "are 0 to 4 (passes of the 1-2-1 filter).");
        }
        // Without surface tension there is no stress to smooth (and no extra ghost cells)
        if (!(sigma_max > 0.0)){
            smoothing = 0;
        }

        // Unsupported combinations
        int coord_sys = 0;
        amrex::ParmParse pp_geom("geometry");
        pp_geom.query("coord_sys", coord_sys);
        if (coord_sys != 0){
            amrex::Abort("SURFACE_TENSION = true is implemented for Cartesian geometry only (geometry.coord_sys = 0): "
                         "the axisymmetric capillary hoop term T_thth/r is missing.");
        }
    };
};

#endif
#endif
