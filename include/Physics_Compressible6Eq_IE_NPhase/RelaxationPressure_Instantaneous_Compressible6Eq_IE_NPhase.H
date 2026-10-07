// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_
#define RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// Iteration cap of the Newton iteration below: it stops after RelaxNewton_MaxIter + 1 = 51
// updates. The test cases need at most 7 iterations; the cap was 100000, which the iteration
// reached where it stalled on its absolute test |f| <= 1e-16 or where the relaxation function has
// no admissible root.
constexpr int RelaxNewton_MaxIter = 50;

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void EvalFx(amrex::Real& fx,
			amrex::Real const& Pe,
			amrex::Array<amrex::Real, NPHASE> const& Alpha,
			amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps,
			Parm const& parm)
{
	fx = 0.0;
	amrex::Real a,b;

	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
		a = 1. + 1./parm.A[i_phase];
		b = parm.B[i_phase]/(parm.A[i_phase]);

		fx += (Alpha_Ps[i_phase] - Alpha[i_phase]*Pe)/(a*Pe +b + 1.e-300);
	}

}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void EvaldFdx(amrex::Real& dfdx,
			amrex::Real const& Pe,
			amrex::Array<amrex::Real, NPHASE> const& Alpha,
			amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps,
			Parm const& parm)
{
	dfdx = 0.0;
	amrex::Real a, b;
	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
		a = 1. + 1./parm.A[i_phase];
		b = parm.B[i_phase]/parm.A[i_phase];

		dfdx += (-Alpha[i_phase]*b - Alpha_Ps[i_phase]*a)/(POW2(a*Pe+b + 1.e-300));
	}

}

// Newton iteration on f(P) = sum_k (Alpha_Ps_k - Alpha_k P)/(gamma_k (P + pinf_k)) (Ndanou et al.
// 2015), from |sum_k Alpha_Ps_k| halved until f >= 0, until |f| <= 1e-16, |f'| < 1e-14 or
// RelaxNewton_MaxIter + 1 updates. Returns the number of updates; the caller checks the result.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
NewtonRaphson(amrex::Real & Pe,
				amrex::Array<amrex::Real, NPHASE> const& Alpha,
				amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps,
				Parm const& parm)
{

	amrex::Real P0 = 0.0;

	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
		P0 += Alpha_Ps[i_phase];
	}

	// searching for initial guess (once P0 has underflowed to 0, halving changes nothing)

	P0 = std::abs(P0);
	amrex::Real f0;
	EvalFx(f0, P0, Alpha, Alpha_Ps, parm);
	int nit = 0;

	while( f0 < 0 && P0 > 0){
		P0 = 0.5*P0;
		EvalFx(f0, P0, Alpha, Alpha_Ps, parm);
		nit++;
		if (nit > 10000){
			break;
		}
	}
	Pe = P0;


	amrex::Real fx, dfdx, eps, epsprime;
	eps = 1.e-16;
	epsprime = 1.e-14;

	EvalFx(fx, Pe, Alpha, Alpha_Ps, parm);

	int niter = 0;
	amrex::Real Pe2;
	while(std::abs(fx) > eps ){
		EvaldFdx(dfdx, Pe, Alpha, Alpha_Ps, parm);
		EvalFx(fx, Pe, Alpha, Alpha_Ps, parm);

		if(std::abs(dfdx) < epsprime){
			break;
		}
		Pe2 = Pe - fx/dfdx;

		Pe = Pe2;
		niter++;

		if(niter > RelaxNewton_MaxIter){
			break;
		}
	}
	return niter;

}

// NASG branch: EOS = 3 with a covolume b_k > 0 or a heat of formation q_k != 0 for some phase.
// EOS = 3 with b = q = 0 takes the stiffened-gas branch, so it reproduces EOS = 1.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool PressureRelaxation_IsNASG (Parm const& parm)
{
	bool NASG = false;
	if (parm.EOS == 3){
		for (int i_phase = 0; i_phase < NPHASE; ++i_phase){
			NASG = NASG || (parm.C[i_phase] > 0.0) || (parm.D[i_phase] != 0.0);
		}
	}
	return NASG;
}

// Instantaneous pressure relaxation of one cell. The volume fractions come from the Newton
// iteration above wherever it does not fail: |f| <= 1e-12 at a finite pressure above every
// -pinf_k of the present phases (alpha_k > 0; the root of the admissible branch), and the sum of
// the volume fractions equal to the sum before to 1e-12 (f carries a 1e-300 guard that the volume
// fractions do not, so |f| alone misses a root at a vanishing P + pinf_k, as for a gas phase at
// zero pressure; the volume fractions also lose digits next to an asymptote). A volume fraction
// out of [0, 1] is not a failure: no volume fraction is clipped, and the bounds come from the
// transport. Where the iteration fails (no convergence, as where a phase energy below its EOS
// minimum leaves no root above the asymptotes, a root beyond an asymptote, a sum not kept, or a
// phase with 0 < alpha_k <= Physics.pressure_relaxation_trace) the volume fractions come from the
// safeguarded relaxation of RelaxationPressure_Safeguarded.H. In both cases
// the phase internal energies are then reset at the pressure given by rho E. Only the volume
// fractions and the phase internal energies change. Iter and Path report the iterations and the
// path (Relax_Path*).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void PressureRelaxation_Cell (amrex::Array<amrex::Real, NSTATE>& ConservativeR,
                              amrex::Array<amrex::Real, NSTATE> const& Conservative,
                              Parm const& parm,
                              bool NASG,
                              int& Iter,
                              int& Path)
{
	amrex::Array<amrex::Real, NSTATE> Primitive;
	Conservative2Primitive(Primitive, Conservative, parm);

	amrex::Array<amrex::Real, NPHASE> Alpha;
	amrex::Array<amrex::Real, NPHASE> Alpha_rho;
	amrex::Array<amrex::Real, NPHASE> Alpha_Ps;
	amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;
	AMREX_D_TERM( amrex::Real u;,
	              amrex::Real v;,
	              amrex::Real w);

	AMREX_D_TERM( amrex::Real rho_u;,
	              amrex::Real rho_v;,
	              amrex::Real rho_w);

	amrex::Real rho, rhoE, rhoe_h;
	amrex::Real P;

	Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), P, Primitive);
	Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u, rho_v, rho_w), rhoE, Conservative);
	amrex::ignore_unused(P, AMREX_D_DECL(rho_u, rho_v, rho_w));

	rho = 0.;
	for(int i_phase=0; i_phase<NPHASE; ++i_phase){
		rho += Alpha_rho[i_phase];
	}

	amrex::Real Pe, gamma, pinf;
	amrex::Array<amrex::Real, NPHASE> AlphaR, Alpha_rho_eR;

	// Partial form of the EOS: phi_k and Pi_k = phi_k p_k. With the NASG EOS (Le Metayer and
	// Saurel 2016), phi_k = alpha_k - b_k alpha_k rho_k and
	// Pi_k = (alpha_k rho_k e_k - alpha_k rho_k q_k - phi_k B_k)/A_k; the phase energies relaxed
	// with the final pressure P give alpha_kR = alpha_k + (Pi_k - phi_k P)/(gamma_k (P + pinf_k)),
	// the stiffened-gas function with alpha_k -> phi_k and alpha_k p_k -> Pi_k. phi_k is taken
	// without the covolume floor, so that a phase compressed past 1/b_k in a mixture cell by the
	// hyperbolic step (phi_k < 0) is relaxed back. Otherwise phi_k = alpha_k and Pi_k = alpha_k p_k.
	amrex::Array<amrex::Real, NPHASE> Phi, Pi;
	for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
		if (NASG){
			Phi[i_phase] = Alpha[i_phase] - parm.C[i_phase]*Alpha_rho[i_phase];
			Pi[i_phase]  = (Alpha_rho_e[i_phase] - Alpha_rho[i_phase]*parm.D[i_phase] - Phi[i_phase]*parm.B[i_phase])/parm.A[i_phase];
		}
		else {
			Phi[i_phase] = Alpha[i_phase];
			Pi[i_phase]  = Alpha_Ps[i_phase];
		}
	}

	int niter = NewtonRaphson(Pe, Phi, Pi, parm);

	for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
		gamma = 1. + 1./parm.A[i_phase];
		pinf = (gamma-1.)*parm.B[i_phase]/gamma;
		if (NASG){
			AlphaR[i_phase] = Alpha[i_phase] + (Pi[i_phase] - Phi[i_phase]*Pe)/(gamma*(pinf+Pe));
		}
		else {
			AlphaR[i_phase] = (Alpha_Ps[i_phase] + Alpha[i_phase]*gamma*pinf + Alpha[i_phase]*Pe*(gamma-1.))/(gamma*(pinf+Pe));
		}
	}

	// Admissibility of the Newton result
	amrex::Real fR;
	EvalFx(fR, Pe, Phi, Pi, parm);
	bool admissible = Relax_Finite(Pe) && (std::abs(fR) <= 1.e-12);
	amrex::Real Sum = 0.0, SumR = 0.0;
	for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
		Sum  += Alpha[i_phase];
		SumR += AlphaR[i_phase];
	}
	admissible = admissible && (std::abs(SumR - Sum) <= 1.e-12);
	for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
		admissible = admissible && Relax_Finite(AlphaR[i_phase]);
		if (Alpha[i_phase] > 0.0){
			const amrex::Real pinf_k = parm.B[i_phase]/(parm.A[i_phase] + 1.0);
			admissible = admissible && (Pe + pinf_k > 0.0)
			          && (Alpha[i_phase] > parm.Physics_Parm.pressure_relaxation_trace);
		}
	}
	if (admissible){
		Iter = niter;
		Path = Relax_PathModel;
	}
	else {
		amrex::Real PR;
		Path = RelaxPressure_Safeguarded<NPHASE>(AlphaR, PR, Iter, Alpha, Alpha_rho, Phi, Pi, parm.A, parm.B,
		                                         parm.Physics_Parm.pressure_relaxation_trace);
	}

	rhoe_h = rhoE - (0.5)*rho*( AMREX_D_TERM( u*u, +v*v, +w*w ) );

	if (NASG){
		// The pressure from the mixture internal energy, with the q_k terms
		amrex::Real sum_BPhi = 0.0, sum_APhi = 0.0, sum_q = 0.0;
		amrex::Array<amrex::Real, NPHASE> PhiR;
		for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
			PhiR[i_phase] = NASG_Xi(AlphaR[i_phase], Alpha_rho[i_phase], parm.C[i_phase])*AlphaR[i_phase];
			sum_BPhi += PhiR[i_phase]*parm.B[i_phase];
			sum_APhi += PhiR[i_phase]*parm.A[i_phase];
			sum_q    += Alpha_rho[i_phase]*parm.D[i_phase];
		}
		amrex::Real pr = (rhoe_h - sum_q - sum_BPhi)/sum_APhi;
		for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
			Alpha_rho_eR[i_phase] = PhiR[i_phase]*(parm.A[i_phase]*pr + parm.B[i_phase]) + Alpha_rho[i_phase]*parm.D[i_phase];
		}
	}
	else {
		// The phase energies are reset at the common pressure pr from
		//   alpha_k rho_k e_k = alpha_k (A_k pr + B_k) + alpha_k rho_k D_k,
		// with alpha_k rho_k unchanged. The D terms (zero for EOS 0 and 1) are added only
		// where D_k != 0, so those EOS keep the same arithmetic.
		amrex::Real pr, a, b, d;
		a = 0.; b = 0.; d = 0.;
		for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
			a += AlphaR[i_phase]*parm.B[i_phase];
			b += AlphaR[i_phase]*parm.A[i_phase];
			if (parm.D[i_phase] != 0.0){
				d += Alpha_rho[i_phase]*parm.D[i_phase];
			}
		}

		pr = (rhoe_h - a - d)/b;
		for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
			Alpha_rho_eR[i_phase] = AlphaR[i_phase]*(parm.A[i_phase]*pr + parm.B[i_phase]);
			if (parm.D[i_phase] != 0.0){
				Alpha_rho_eR[i_phase] += Alpha_rho[i_phase]*parm.D[i_phase];
			}
		}
	}

	ConservativeR = Conservative;
	for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
		ConservativeR[INDEX_VolumeFraction1 + i_phase] = AlphaR[i_phase];
		ConservativeR[INDEX_Energy1 + i_phase]         = Alpha_rho_eR[i_phase];
	}
}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
PressureRelaxationInstantaneous (amrex::MultiFab & mf,
							     const Parm parm)
{
	const bool NASG = PressureRelaxation_IsNASG(parm);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {


            // Get the box
            const Box& bx = mfi.tilebox();

            // Pull the data into an array
            Array4<Real const> statein   = mf.const_array(mfi);
            Array4<Real      > stateout  = mf.array(mfi);


            // ===========================
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
            	amrex::Array<amrex::Real, NSTATE> Conservative;
            	amrex::Array<amrex::Real, NSTATE> ConservativeR;

            	for (int v = 0; v < NSTATE; v++){
            		Conservative[v] = statein(i,j,k,v);
            	}

				int Iter, Path;
				PressureRelaxation_Cell(ConservativeR, Conservative, parm, NASG, Iter, Path);

				for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
					stateout(i,j,k, INDEX_VolumeFraction1 + i_phase) = ConservativeR[INDEX_VolumeFraction1 + i_phase];
					stateout(i,j,k, INDEX_Energy1 + i_phase)         = ConservativeR[INDEX_Energy1 + i_phase];
				}

			});
		}
	}

}

#endif
#endif
