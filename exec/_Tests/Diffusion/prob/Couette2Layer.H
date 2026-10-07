// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Two-layer Couette flow
// A problem of the Diffusion case (prob.problem = Couette2Layer), included by ProblemICBC.H,
// for the models 5Eq, 6Eq and 6Eq-N (NPHASE >= 2, phases 3 and up absent).

#ifndef DIFFUSION_COUETTE2LAYER_H_
#define DIFFUSION_COUETTE2LAYER_H_

#ifdef DIFFUSION_HAS_COUETTE2LAYER
namespace Couette2Layer {

// =======================================
// Initial condition
//
// Two-layer Couette flow: phase 1 (viscosity mu_1) fills y_lo < y < h and phase 2 (mu_2)
// fills h < y < y_hi, the wall y = y_lo is at rest and the wall y = y_hi moves at the speed U
// along x. The steady flow has the uniform shear stress tau = U/((h - y_lo)/mu_1 + (y_hi - h)/mu_2)
// and the piecewise linear velocity
//   u = tau (y - y_lo)/mu_1                    for y < h,
//   u = tau (h - y_lo)/mu_1 + tau (y - h)/mu_2 for y > h,
// the initial condition. The pressure is the uniform prob.P0 and the phase densities are
// prob.rho_1 and prob.rho_2 (prob.rho in the N-phase model, where phases 3 and up are absent);
// with U well below the sound speeds the flow stays incompressible.
// =======================================

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real CouetteVelocity(amrex::Real y, Parm const& parm)
{
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    return (y < parm.h) ? parm.tau*(y - parm.y_lo)/parm.mu_1
                        : parm.tau*(parm.h - parm.y_lo)/parm.mu_1 + parm.tau*(y - parm.h)/parm.mu_2;
#else
    return (y < parm.h) ? parm.tau*(y - parm.y_lo)/parm.mu[0]
                        : parm.tau*(parm.h - parm.y_lo)/parm.mu[0] + parm.tau*(y - parm.h)/parm.mu[1];
#endif
}

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
    amrex::ignore_unused(x, z, AMREX_D_DECL(dx, dy, dz),
                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));

#if (AMREX_SPACEDIM > 1)
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier*dy_FinestLevel;
    const amrex::Real H   = (eta > 0.0) ? 0.5*(1.0 - std::tanh((y - parm.h)/(std::sqrt(2.0)*eta)))
                                        : ((y < parm.h) ? 1.0 : 0.0);
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    const amrex::Real Alpha1 = parm.eps + (1.0 - 2.0*parm.eps)*H;
    const amrex::Real Alpha2 = 1.0 - Alpha1;
#endif

    const amrex::Real u = CouetteVelocity(y, parm);
    const amrex::Real v = 0.0;
    const amrex::Real w = 0.0;
    amrex::ignore_unused(w); // used in 3D only

#if (PHYSICS == FIVEEQS)
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*parm.rho_1, Alpha2*parm.rho_2,
       AMREX_D_DECL(u, v, w),
                    parm.P0,
                    Alpha1);
#elif (PHYSICS == SIXEQS)
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha1*parm.rho_1, Alpha2*parm.rho_2,
       AMREX_D_DECL(u, v, w),
                    Alpha1*parm.P0, Alpha2*parm.P0,
                    Alpha1);
#else
    const amrex::Real P = parm.P0;
    amrex::Array<amrex::Real, NPHASE> Alpha, Alpha_rho, Alpha_P;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha[i_Phase] = 0.0;
    }
    Alpha[0] = parm.eps + (1.0 - 2.0*parm.eps)*H;
    Alpha[1] = 1.0 - Alpha[0];
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho[i_Phase] = Alpha[i_Phase]*parm.rho_phase[i_Phase];
        Alpha_P  [i_Phase] = Alpha[i_Phase]*P;
    }
    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alpha, Alpha_rho, Alpha_P,
       AMREX_D_DECL(u, v, w),
                    P);
#endif
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
#else
    amrex::ignore_unused(y, Conservative, parm);
    amrex::Abort("Couette2Layer needs DIM = 2 or 3.");
#endif
}


// =======================================
// Boundary condition (run.lo_bc and run.hi_bc = 1 in y)
//
// No-slip walls: the wall y = y_lo is at rest and the wall y = y_hi moves at (U, 0, 0). The
// ghost cell takes the state of its mirror cell across the wall with the velocity reflected
// about the wall velocity, u_g = 2 u_wall - u, and the kinetic energy in rho E (in each phase
// energy in the six-equation model) replaced by that of u_g, so the face velocity is u_wall and
// the face pressure, volume fractions, densities and phase internal energies are those of the
// mirror cell. Other boundaries: zero gradient.
// =======================================

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
    amrex::ignore_unused(xBC, yBC, zBC, time, AMREX_D_DECL(ilo, jlo, klo), AMREX_D_DECL(ihi, jhi, khi),
                         AMREX_D_DECL(dx, dy, dz));

#if (AMREX_SPACEDIM > 1)
    if (j < jlo || j > jhi)
    {
        const int  jm     = (j < jlo) ? 2*jlo - j - 1 : 2*jhi - j + 1;
        const Real U_wall = (j < jlo) ? 0.0 : parm->U;

        Array<Real,NSTATE> Interior;
        for (int iState = 0; iState < NSTATE; ++iState){
            Interior[iState] = q(i,jm,k,iState);
        }
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
        const Real rho = Interior[INDEX_Mass1] + Interior[INDEX_Mass2];
#else
        Real rho = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            rho += Interior[INDEX_Mass1 + i_Phase];
        }
#endif
        AMREX_D_TERM( const Real u = Interior[INDEX_MomentumX]/rho;,
                      const Real v = Interior[INDEX_MomentumY]/rho;,
                      const Real w = Interior[INDEX_MomentumZ]/rho );
        AMREX_D_TERM( const Real u_g = 2.0*U_wall - u;,
                      const Real v_g = -v;,
                      const Real w_g = -w );
        const Real dKE = 0.5*(AMREX_D_TERM( u_g*u_g - u*u,
                                          + v_g*v_g - v*v,
                                          + w_g*w_g - w*w ));

        Array<Real,NSTATE> Ghost = Interior;
        AMREX_D_TERM( Ghost[INDEX_MomentumX] = rho*u_g;,
                      Ghost[INDEX_MomentumY] = rho*v_g;,
                      Ghost[INDEX_MomentumZ] = rho*w_g );
#if (PHYSICS == FIVEEQS)
        Ghost[INDEX_Energy] += rho*dKE;
#elif (PHYSICS == SIXEQS)
        Ghost[INDEX_Energy1] += Interior[INDEX_Mass1]*dKE;
        Ghost[INDEX_Energy2] += Interior[INDEX_Mass2]*dKE;
#else
        Ghost[INDEX_Total_Energy] += rho*dKE;
#endif

        // Dirichlet on the face value, (ghost + mirror)/2
        BC_Dirichlet = 1.0;
        BC_Neumann   = 0.0;
        BC_Function  = 0.5*(Ghost[n] + Interior[n]);
        return;
    }
#else
    amrex::ignore_unused(i, j, k, n, q, parm);
#endif
    // Zero gradient
    BC_Dirichlet = 0.0;
    BC_Neumann   = 1.0;
    BC_Function  = 0.0;
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/couette.csv with, over the composite grid,
//   errU_L1, errU_max   int |u - u_e| dV/(V U) and max |u - u_e|/U, with u_e the exact
//                       steady profile at the cell center and V the area of the box
//   tau_lo, tau_hi      the wall shear stress at y_lo and y_hi, the mean over x of the
//                       stress of the wall face, mu (u_1 - u_wall)/(dy/2) up to its sign, with
//                       u_1 the velocity of the wall cell and mu its mixture viscosity
//   tau_e               the exact shear stress
//   mass_1, mass_2      int alpha_k rho_k dV
//   energy              int rho E dV (the sum of the phase energies in the six-equation model)
//   P_dev               max |P - P0|/P0
// =======================================

template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const auto v = Measure(o, {Sum, Max, Sum, Sum, Sum, Sum, Sum, Max, Sum}, 0,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const& parm) -> Values<9>
    {
        amrex::Array<Real, NSTATE> Cons;
        for (int n = 0; n < NSTATE; ++n) { Cons[n] = U(c.i,c.j,c.k,n); }
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
        const Real rho = Cons[INDEX_Mass1] + Cons[INDEX_Mass2];
        const Real u   = Cons[INDEX_MomentumX]/rho;
#else
        Real rho = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            rho += Cons[INDEX_Mass1 + i_Phase];
        }
        const Real u = Cons[INDEX_MomentumX]/rho;
#endif
        Real P;
        Conservative2Pressure(P, Cons, parm);
        Real mu, mub;
        Conservative2Viscosity(mu, mub, Cons, parm);
#if (AMREX_SPACEDIM > 1)
        const Real err    = std::abs(u - CouetteVelocity(c.x[1], parm))/parm.U;
        const Real wall   = c.dx[0]/parm.L_x/(0.5*c.dx[1]);
        const Real tau_lo = (c.j == c.domain.smallEnd(1)) ? mu*u*wall : 0.0;
        const Real tau_hi = (c.j == c.domain.bigEnd(1))   ? mu*(parm.U - u)*wall : 0.0;
#else
        const Real err = 0.0, tau_lo = 0.0, tau_hi = 0.0;
#endif
#if (PHYSICS == FIVEEQS)
        return {err*c.dV, err, tau_lo, tau_hi, Cons[INDEX_Mass1]*c.dV, Cons[INDEX_Mass2]*c.dV,
                Cons[INDEX_Energy]*c.dV, std::abs(P - parm.P0)/parm.P0, c.dV};
#elif (PHYSICS == SIXEQS)
        return {err*c.dV, err, tau_lo, tau_hi, Cons[INDEX_Mass1]*c.dV, Cons[INDEX_Mass2]*c.dV,
                (Cons[INDEX_Energy1] + Cons[INDEX_Energy2])*c.dV, std::abs(P - parm.P0)/parm.P0, c.dV};
#else
        return {err*c.dV, err, tau_lo, tau_hi, Cons[INDEX_Mass1]*c.dV, Cons[INDEX_Mass1 + 1]*c.dV,
                Cons[INDEX_Total_Energy]*c.dV, std::abs(P - parm.P0)/parm.P0, c.dV};
#endif
    });
    Parm const& hp = HostParm(o);
    WriteRow(o, "couette.csv", {{"errU_L1", v[0]/v[8]}, {"errU_max", v[1]},
                                {"tau_lo", v[2]}, {"tau_hi", v[3]}, {"tau_e", hp.tau},
                                {"mass_1", v[4]}, {"mass_2", v[5]}, {"energy", v[6]},
                                {"P_dev", v[7]}});
}

} // namespace Couette2Layer
#endif // DIFFUSION_HAS_COUETTE2LAYER

#endif // DIFFUSION_COUETTE2LAYER_H_
