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

				NewtonRaphson(Pe, Alpha, Alpha_Ps, parm);
			


                amrex::Real suma = 0.;

                for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
                    gamma = 1. + 1./parm.A[i_phase];
                    pinf = (gamma-1.)*parm.B[i_phase]/gamma;
                                        
                    AlphaR[i_phase] = (Alpha_Ps[i_phase] + Alpha[i_phase]*gamma*pinf + Alpha[i_phase]*Pe*(gamma-1.))/(gamma*(pinf+Pe));
                    
                    suma += AlphaR[i_phase];
                }

				
				rhoe_h = rhoE - (0.5)*rho*( AMREX_D_TERM( u*u, +v*v, +w*w ) );
		
				
				amrex::Real pr, a, b;
				a = 0.; b = 0.;
				for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
					a += AlphaR[i_phase]*parm.B[i_phase];
					b += AlphaR[i_phase]*parm.A[i_phase];
				}
				

				pr = (rhoe_h - a)/b;
				for(int i_phase = 0 ; i_phase<NPHASE; ++i_phase){
					Alpha_rho_eR[i_phase] = AlphaR[i_phase]*(parm.A[i_phase]*pr + parm.B[i_phase]); 
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
