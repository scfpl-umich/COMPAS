// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <AMReX_Gpu.H>
#include <Macros.H>
#include <LoopMacros.H>
#include <Parm.H>

#include <string>

namespace {

// Names the value that failed a check. i_Phase >= 0 is the entry of an N-phase list.
std::string EOSValueName (std::string const& key, int i_Phase)
{
    return (i_Phase < 0) ? key : "value " + std::to_string(i_Phase+1) + " of " + key;
}

// An EOS key that the inputs file does not set is left at 0, which gamma > 1 rejects.
void CheckGamma (amrex::Real gamma, std::string const& key, int EOS, int i_Phase = -1)
{
    if (!(gamma > 1.0)) {
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " requires gamma > 1, but " + EOSValueName(key, i_Phase)
                     + " is " + std::to_string(gamma) + ". Set " + key + " in the inputs file.");
    }
}

// The temperature is (rho e - rho D - B_T)/(rho cv), so it needs cv > 0.
[[maybe_unused]] void CheckCv (amrex::Real cv, std::string const& key, std::string const& user, int i_Phase = -1)
{
    if (!(cv > 0.0)) {
        amrex::Abort(user + " uses the temperature, which requires cv > 0, but " + EOSValueName(key, i_Phase)
                     + " is " + std::to_string(cv) + ". Set " + key + " in the inputs file.");
    }
}

#if (PHYSICS == FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)
// Each per-phase EOS list needs exactly NPHASE values. A shorter list leaves the last phases unset.
void CheckPhaseLists ()
{
    amrex::ParmParse pp("prob.EOS");
    for (std::string const key : {"gamma", "pinf", "gGamma", "pref", "eref", "cv", "mu", "muB", "kappa"}){
        if (pp.contains(key.c_str()) && pp.countval(key.c_str()) != NPHASE){
            amrex::Abort("prob.EOS." + key + " has " + std::to_string(pp.countval(key.c_str()))
                         + " values, but NPHASE = " + std::to_string(NPHASE) + " needs one value per phase.");
        }
    }
}
#endif

}



#if (PHYSICS == FIVEEQS)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    if (EOS == 0 || EOS == 1){
        CheckGamma(gamma_1, "prob.EOS.gamma_1", EOS);
        CheckGamma(gamma_2, "prob.EOS.gamma_2", EOS);
    }
#if (DIFFUSION == true)
    CheckCv(cv_1, "prob.EOS.cv_1", "DIFFUSION = true");
    CheckCv(cv_2, "prob.EOS.cv_2", "DIFFUSION = true");
#endif

    switch(EOS){
    case 0:  // Ideal gas
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = 0.; B2 = 0.;
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 1: // SG
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 2: // Mie-Gruneisen
        if (gGamma_1 <= 0.0 || gGamma_2 <= 0.0){
            amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma_1 > 0 and prob.EOS.gGamma_2 > 0.");
        }
        A1 = 1./gGamma_1; A2 = 1./gGamma_2;
        B1 = -pref_1/gGamma_1; B2 = -pref_2/gGamma_2;
        C1 = 0.; C2 = 0.; 
        D1 = eref_1; D2 = eref_2;
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen).");
        break;
    }

    if (Physics_Parm.source_term == 0){
        p_min[0] = -B1/A1;
        p_min[1] = -B2/A2;
    }
    else if (Physics_Parm.source_term == 1){
        p_min[0] = -B1/(A1 + 1.0);
        p_min[1] = -B2/(A2 + 1.0);
    }
    else {
        amrex::Abort("Physics.source_term must be 0 or 1.");
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {inf, inf, AMREX_D_DECL( inf,  inf,  inf),   inf, 1.0};
    LowerBound = {0.0, 0.0, AMREX_D_DECL(-inf, -inf, -inf),  -inf, 0.0};
    ID_BoundType   = { -1,  -1, AMREX_D_DECL(   2,    2,    2),    2,   0};

}

#endif

#if (PHYSICS == FIVEEQS_NPHASE)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    CheckPhaseLists();
    if (EOS == 0 || EOS == 1){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            CheckGamma(gamma[i_Phase], "prob.EOS.gamma", EOS, i_Phase);
        }
    }
#if (DIFFUSION == true)
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        CheckCv(cv[i_Phase], "prob.EOS.cv", "DIFFUSION = true", i_Phase);
    }
#endif

    switch(EOS){
    case 0:  // Ideal gas
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = 0.0;
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 1: // SG
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 2: // Mie-Gruneisen
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (gGamma[i_Phase] <= 0.0){
                amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma > 0 for every phase.");
            }
            A[i_Phase] = 1./gGamma[i_Phase];
            B[i_Phase] = -pref[i_Phase]/gGamma[i_Phase];
            C[i_Phase] = 0.0;
            D[i_Phase] = eref[i_Phase];
        }
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen).");
        break;
    }

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        if (Physics_Parm.source_term == 0){
            p_min[i_Phase] = -B[i_Phase]/A[i_Phase];
        }
        else if (Physics_Parm.source_term == 1){
            p_min[i_Phase] = -B[i_Phase]/(A[i_Phase] + 1.0);
        }
        else {
            amrex::Abort("Physics.source_term must be 0 or 1.");
        }
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {D_LOOP(NPHASE,inf), AMREX_D_DECL( inf,  inf,  inf),   inf, D_LOOP(NPHASE,1.0)};
    LowerBound = {D_LOOP(NPHASE,0.0), AMREX_D_DECL(-inf, -inf, -inf),  -inf, D_LOOP(NPHASE,0.0)};
    ID_BoundType   = {D_LOOP(NPHASE,-1),  AMREX_D_DECL(   2,    2,    2),    2,   D_LOOP(NPHASE,0)};

}

#endif


#if (PHYSICS == SIXEQS)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    if (EOS == 0 || EOS == 1){
        CheckGamma(gamma_1, "prob.EOS.gamma_1", EOS);
        CheckGamma(gamma_2, "prob.EOS.gamma_2", EOS);
    }
    if (Physics_Parm.pressure_temperature_relaxation){
        CheckCv(cv_1, "prob.EOS.cv_1", "Physics.pressure_temperature_relaxation");
        CheckCv(cv_2, "prob.EOS.cv_2", "Physics.pressure_temperature_relaxation");
    }

    switch(EOS){
    case 0:  // Ideal gas
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = 0.; B2 = 0.;
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 1: // SG
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
     case 2: // Mie-Gruneisen
        if (gGamma_1 <= 0.0 || gGamma_2 <= 0.0){
            amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma_1 > 0 and prob.EOS.gGamma_2 > 0.");
        }
        A1 = 1./gGamma_1; A2 = 1./gGamma_2;
        B1 = -pref_1/gGamma_1; B2 = -pref_2/gGamma_2;
        C1 = 0.; C2 = 0.; 
        D1 = eref_1; D2 = eref_2;
        break;
     default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen).");
        break;
    }


    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    p_min[0] = -B1/(A1 + 1.0);
    p_min[1] = -B2/(A2 + 1.0);
    UpperBound = {inf, inf, AMREX_D_DECL( inf,  inf,  inf),     inf,     inf,  1.0};
    LowerBound = {0.0, 0.0, AMREX_D_DECL(-inf, -inf, -inf),    -inf,    -inf,  0.0};
    ID_BoundType   = { -1,  -1, AMREX_D_DECL(   2,    2,    2),      2,      2,    0};

}

#endif



#if (PHYSICS == SIXEQS_IE_NPHASE)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    CheckPhaseLists();
    if (EOS == 0 || EOS == 1){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            CheckGamma(gamma[i_Phase], "prob.EOS.gamma", EOS, i_Phase);
        }
    }
#if (DIFFUSION == true)
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        CheckCv(cv[i_Phase], "prob.EOS.cv", "DIFFUSION = true", i_Phase);
    }
#endif

    switch(EOS){
    case 0:  // Ideal gas
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = 0.0;
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 1: // SG
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 2: // Mie-Gruneisen
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (gGamma[i_Phase] <= 0.0){
                amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma > 0 for every phase.");
            }
            A[i_Phase] = 1./gGamma[i_Phase];
            B[i_Phase] = -pref[i_Phase]/gGamma[i_Phase];
            C[i_Phase] = 0.0;
            D[i_Phase] = eref[i_Phase];
        }
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen).");
        break;
    }

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
       
        p_min[i_Phase] = -B[i_Phase]/(A[i_Phase] + 1.0);
    
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {D_LOOP(NPHASE,1.0), D_LOOP(NPHASE,inf), D_LOOP(NPHASE,inf), AMREX_D_DECL( inf,  inf,  inf),   inf };
    LowerBound = {D_LOOP(NPHASE,0.0), D_LOOP(NPHASE,0.0), D_LOOP(NPHASE,0.0), AMREX_D_DECL(-inf, -inf, -inf),  -inf};
    ID_BoundType   = {D_LOOP(NPHASE,0), D_LOOP(NPHASE,-1), D_LOOP(NPHASE,-1),  AMREX_D_DECL(   2,    2,    2),    2};

}

#endif