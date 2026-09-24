// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>


// =======================================
// Initial condition
//
// The box is filled with fluid 1 (alpha_1 = 1) at rest at pressure P0. Fluid 2 enters
// through the inflow boundary set in BoundaryFunction_Robin.
// =======================================

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
{
    // Prevent unused-variable warnings
    sharp = false;
    amrex::ignore_unused(x, y, z, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

    // ---------------------------------------
    // Fluid 1 at rest
    // ---------------------------------------
    const amrex::Real Alpha1 = 1.0;
    const amrex::Real Alpha2 = 1.0 - Alpha1;

    const amrex::Real Alpha1_rho1 = Alpha1 * parm.rho1;
    const amrex::Real Alpha2_rho2 = Alpha2 * parm.rho2;

    const amrex::Real u = 0.0;
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    const amrex::Real P = parm.P0;

    // ---------------------------------------
    // Convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1_rho1, Alpha2_rho2,
                    AMREX_D_DECL(u, v, w),
                    P,
                    Alpha1);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// User boundary condition, in Robin form
//
// COMPAS calls this for the ghost cells of each boundary of type 1 (Inflow), here x = 0,
// and fills each conservative component n from
//   BC_Dirichlet q_bc + BC_Neumann dq/dn = BC_Function,
// with q_bc the mean of the ghost cell and its mirror interior cell. At x = 0, where
// r = sqrt(y^2 + z^2) < r_jet, fluid 2 enters: the boundary value of the partial densities
// and volume fraction is the mean of the first interior cell and pure fluid 2, the velocity
// is (u_jet, 0, 0) and the pressure is that of the first interior cell. Elsewhere on x = 0
// the boundary is a slip wall (u = 0). Other faces get a zero-gradient condition, used
// only if their type is set to 1.
// =======================================

#define USER_BOUNDARY_FUNC_ROBIN
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_Robin(int i, int j, int k, int n,
                            Real& BC_Dirichlet,
                            Real& BC_Neumann,
                            Real& BC_Function,
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
{
    amrex::ignore_unused(xBC, time, AMREX_D_DECL(ilo, jlo, klo), AMREX_D_DECL(ihi, jhi, khi),
                         AMREX_D_DECL(dx, dy, dz));

    constexpr Real r_jet = 0.1;   // radius of the inflow opening

    if (i < ilo)
    {
        // State of the first interior cell
        Array<Real,NSTATE> Conservative, Primitive;
        for (int iState = 0; iState < NSTATE; ++iState){
            Conservative[iState] = q(ilo,j,k,iState);
        }
        Real Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u, v, w), p, Alpha1;
        Conservative2Primitive(Primitive, Conservative, *parm);
        Primitive2State(Alpha1_rho1, Alpha2_rho2, AMREX_D_DECL(u, v, w), p, Alpha1, Primitive);

        // Boundary state
        Real Alpha1_rho1_bc, Alpha2_rho2_bc, AMREX_D_DECL(u_bc, v_bc, w_bc), Alpha1_bc;
        const Real rBC = std::sqrt(zBC*zBC + yBC*yBC);
        if (rBC < r_jet){
            // Inflow of pure fluid 2
            const Real Alpha1_Ghost = 0.0;
            const Real Alpha2_Ghost = 1.0 - Alpha1_Ghost;
            Alpha1_rho1_bc = (Alpha1_rho1 + Alpha1_Ghost*parm->rho1)/2.0;
            Alpha2_rho2_bc = (Alpha2_rho2 + Alpha2_Ghost*parm->rho2)/2.0;
            AMREX_D_TERM(u_bc = parm->u_jet;,
                         v_bc = 0.0;,
                         w_bc = 0.0);
            Alpha1_bc = (Alpha1 + Alpha1_Ghost)/2.0;
        }
        else{
            // Slip wall
            Alpha1_rho1_bc = Alpha1_rho1;
            Alpha2_rho2_bc = Alpha2_rho2;
            AMREX_D_TERM(u_bc = 0.0;,
                         v_bc = v;,
                         w_bc = w);
            Alpha1_bc = Alpha1;
        }
        const Real p_bc = p;

        Array<Real,NSTATE> Conservative_bc, Primitive_bc;
        State2Primitive(Primitive_bc, Alpha1_rho1_bc, Alpha2_rho2_bc, AMREX_D_DECL(u_bc, v_bc, w_bc), p_bc, Alpha1_bc);
        Primitive2Conservative(Conservative_bc, Primitive_bc, *parm);

        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = Conservative_bc[n];
    }
    else
    {
        // Zero gradient
        BC_Dirichlet = 0.0;
        BC_Neumann   = 1.0;
        BC_Function  = 0.0;
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/jet.csv with
//   x_tip   largest x of a cell center with alpha_2 >= 0.5, the penetration of the jet
//   V_2     volume of fluid 2, int alpha_2 dV (an area per unit depth in 2D)
//   m_2     mass of fluid 2, int alpha_2 rho_2 dV
// The maximum and the integrals are over the uncovered cells of the composite grid.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;
    const auto v = Measure(o, {Max, Sum, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<3>
    {
        const Real a2 = 1.0 - U(c.i,c.j,c.k,INDEX_VolumeFraction1);
        return {(a2 >= 0.5) ? c.x[0] : -big, a2 * c.dV, U(c.i,c.j,c.k,INDEX_Mass2) * c.dV};
    });
    WriteRow(o, "jet.csv", {{"x_tip", (v[0] > -big) ? v[0] : 0.0}, {"V_2", v[1]}, {"m_2", v[2]}});
}


#endif
