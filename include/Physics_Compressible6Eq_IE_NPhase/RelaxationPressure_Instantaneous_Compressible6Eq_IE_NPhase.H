// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_
#define RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void EvalFx(amrex::Real& fx,
			amrex::Real const& Pe, 
			amrex::Array<amrex::Real, NPHASE> const& Alpha, 
			amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps, 
			const Parm parm)
{
	fx = 0.0;
	amrex::Real a,b;

	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
	// for (int i_phase : iterL) {
		a = 1. + 1./parm.A[i_phase];
		b = parm.B[i_phase]/(parm.A[i_phase]);

		// fx += (Alpha_Ps[i_phase] - Alpha[i_phase]*Pe)/(std::max(a*Pe +b, 1.e-15));
		fx += (Alpha_Ps[i_phase] - Alpha[i_phase]*Pe)/(a*Pe +b + 1.e-300);
	}

}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void EvaldFdx(amrex::Real& dfdx,
			amrex::Real const& Pe, 
			amrex::Array<amrex::Real, NPHASE> const& Alpha, 
			amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps, 
			const Parm parm)
{
	dfdx = 0.0;
	amrex::Real a, b;
	// for (int i_phase : iterL) {
	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
		a = 1. + 1./parm.A[i_phase];
		b = parm.B[i_phase]/parm.A[i_phase];

		// dfdx += (-Alpha[i_phase]*b - Alpha_Ps[i_phase]*a)/(std::max(POW2(a*Pe+b), 1.e-15));
		dfdx += (-Alpha[i_phase]*b - Alpha_Ps[i_phase]*a)/(POW2(a*Pe+b + 1.e-300));
	}

}

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void 
NewtonRaphson(amrex::Real & Pe, 
				amrex::Array<amrex::Real, NPHASE> const& Alpha, 
				amrex::Array<amrex::Real, NPHASE> const& Alpha_Ps, 
				const Parm parm)
{

	
	// Pe = 0.;
	amrex::Real P0 = 0.0;

	for (int i_phase=0; i_phase< NPHASE; ++i_phase){
		P0 += Alpha_Ps[i_phase];
	}

	// searching for initial guess

	P0 = std::abs(P0);
	amrex::Real f0;
	EvalFx(f0, P0, Alpha, Alpha_Ps, parm);
	int nit = 0;

	while( f0 < 0){
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

	int MaxIter = 100000;
	int niter = 0;
	amrex::Real Pe2;
	while(std::abs(fx) > eps ){
	// while(fx > 0. ){
		EvaldFdx(dfdx, Pe, Alpha, Alpha_Ps, parm);
		EvalFx(fx, Pe, Alpha, Alpha_Ps, parm);

		if(std::abs(dfdx) < epsprime){
			break;
		}
		Pe2 = Pe - fx/dfdx;

		// if(Pe2 < Pmin){
		// 	Pe2 = Pmin + epsP*niter*(Pmax-Pmin);
		// }
		// if(Pe2 > Pmax){
		// 	Pe2 = Pmax - epsP*niter*(Pmax-Pmin);
		// }
		
		Pe = Pe2;
		niter++;

		if(niter > 100000){
			amrex::Print()<<"niter: "<<niter<<" fx: "<<fx<<" dfdx: "<<dfdx<<"\n";
			break;
		}
	}

}


AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
PressureRelaxationInstantaneous (amrex::MultiFab & mf,
							     const Parm parm)
{
	// NASG branch: EOS = 3 with a covolume b_k > 0 or a heat of formation q_k != 0 for some phase.
	// EOS = 3 with b = q = 0 takes the stiffened-gas branch, so it reproduces EOS = 1.
	bool NASG = false;
	if (parm.EOS == 3){
		for (int i_phase = 0; i_phase < NPHASE; ++i_phase){
			NASG = NASG || (parm.C[i_phase] > 0.0) || (parm.D[i_phase] != 0.0);
		}
	}

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {


            // Get the box and the box with a ghost layer
            const Box& bx = mfi.tilebox();

            // Pull the data into an array
            Array4<Real const> statein   = mf.const_array(mfi);
            Array4<Real      > stateout  = mf.array(mfi);

            
            // ===========================
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
            	amrex::Array<amrex::Real, NSTATE> Conservative;
            	amrex::Array<amrex::Real, NSTATE> Primitive;

				// amrex::Print()<<"Here\n";

            	for (int v = 0; v < NSTATE; v++){
            		Conservative[v] = statein(i,j,k,v);
					// amrex::Print()<<"Conservative "<<v<<" "<<Conservative[v]<<"\n";
            	}
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

                amrex::Real rho, rhoE, rhoe_h, rhoe_e;
                amrex::Real P;
				amrex::Real eps = 0.*1.e-10;

                Primitive2State(Alpha, Alpha_rho, Alpha_Ps, AMREX_D_DECL(u,v,w), P, Primitive);
				Conservative2State(Alpha, Alpha_rho, Alpha_rho_e, AMREX_D_DECL(rho_u, rho_v, rho_w), rhoE, Conservative);

				rho = 0.;
				for(int i_phase=0; i_phase<NPHASE; ++i_phase){
					rho += Alpha_rho[i_phase];
				}
				
				amrex::Real Pe, gamma, pinf;
				amrex::Real np = (double) NPHASE;
				amrex::Array<amrex::Real, NPHASE> AlphaR, Alpha_rho_eR;

				if (NASG){
					// NASG (Le Metayer and Saurel 2016). With phi_k = alpha_k - b_k alpha_k rho_k and
					// Pi_k = phi_k p_k = (alpha_k rho_k e_k - alpha_k rho_k q_k - phi_k B_k)/A_k, the
					// phase energies relaxed with the final pressure P give
					//   alpha_kR = alpha_k + (Pi_k - phi_k P)/(gamma_k (P + pinf_k)),
					// and sum_k alpha_kR = 1 is f(P) = sum_k (Pi_k - phi_k P)/(gamma_k P + gamma_k pinf_k) = 0,
					// the stiffened-gas function with alpha_k -> phi_k and alpha_k p_k -> Pi_k. The
					// pressure is then recomputed from the mixture internal energy, with the q_k terms.
					// phi_k is taken without the covolume floor, so that a phase compressed past 1/b_k
					// in a mixture cell by the hyperbolic step (phi_k < 0) is relaxed back.
					amrex::Array<amrex::Real, NPHASE> Phi, Pi;
					for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
						Phi[i_phase] = Alpha[i_phase] - parm.C[i_phase]*Alpha_rho[i_phase];
						Pi[i_phase]  = (Alpha_rho_e[i_phase] - Alpha_rho[i_phase]*parm.D[i_phase] - Phi[i_phase]*parm.B[i_phase])/parm.A[i_phase];
					}
					NewtonRaphson(Pe, Phi, Pi, parm);

					for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
						gamma = 1. + 1./parm.A[i_phase];
						pinf = (gamma-1.)*parm.B[i_phase]/gamma;
						AlphaR[i_phase] = Alpha[i_phase] + (Pi[i_phase] - Phi[i_phase]*Pe)/(gamma*(pinf+Pe));
					}

					rhoe_h = rhoE - (0.5)*rho*( AMREX_D_TERM( u*u, +v*v, +w*w ) );

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
						stateout(i,j,k, i_phase) = AlphaR[i_phase];
						stateout(i,j,k, 2*NPHASE + i_phase) = Alpha_rho_eR[i_phase];
					}
					return;
				}

				NewtonRaphson(Pe, Alpha, Alpha_Ps, parm);
			


                amrex::Real suma = 0.;

                for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
                    gamma = 1. + 1./parm.A[i_phase];
                    pinf = (gamma-1.)*parm.B[i_phase]/gamma;
                                        
                    AlphaR[i_phase] = (Alpha_Ps[i_phase] + Alpha[i_phase]*gamma*pinf + Alpha[i_phase]*Pe*(gamma-1.))/(gamma*(pinf+Pe));
                    
                    suma += AlphaR[i_phase];
                }

				
				rhoe_h = rhoE - (0.5)*rho*( AMREX_D_TERM( u*u, +v*v, +w*w ) );
		
				
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


				for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
					stateout(i,j,k, i_phase) = AlphaR[i_phase];
					stateout(i,j,k, 2*NPHASE + i_phase) = Alpha_rho_eR[i_phase];
				}

			});
		}
	}

}

#endif
#endif
