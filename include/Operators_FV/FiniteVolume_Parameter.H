// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FiniteVolume_Parameter_H_
#define FiniteVolume_Parameter_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>

enum class FVReconstruction : int {
    GODUNOV,
    MUSCL_MM,
    MUSCL_VL,
    MUSCL_SB,
    MUSCL_MC,
    WENO_3,
    WENO_3B,
    WENO_5,
    WENO_5B
};

struct FiniteVolume_Parameter
{
    FVReconstruction Scheme = FVReconstruction::MUSCL_MM;
    int Num_Grow   = 2;
    int ID_WENO    = 0;
    int ID_Limiter = 1;
    int ID_Poly    = 1;
    int ID_Quad    = 0;

    int ID_THINC    = 0;
    int ID_THINC_UseNormal = 1;
    amrex::Real THINC_EPS = 1.0e-4;
    amrex::Real THINC_Beta = 2.3;
    amrex::Real THINC_NormalThreshold = 0.25;
    int THINC_Multi = 0;        //0: THINC along each direction; 1: multi-dimensional THINC (MTHINC)
    int THINC_FaceState = 0;    //0: MUSCL-minmod state on both sides of a THINC face; 1: state of the base scheme, only the THINC side is rebuilt
    int THINC_PhaseDensity = 0; //0: phase densities of a THINC side from its face values; 1: from the cell
    int ID_THINC_Options = 0;   //1 when Multi, FaceState or PhaseDensity is not 0 (set in Initialize)


    int ID_Bound = 0; //0: Bound-preserving is off; 1: Bound-preserving is on
    int ID_Return = 0; //1: restart the step. Set and read on the host only, through h_parm
    int ID_Lambda = 1; //0: General; 1: Special (like all gases or single-phase flows)
    amrex::Real Quadrature_Weight = 0.5; //cmax*dt/dx/Quadrature_Weight=CFL

    //Read parameters from the inputs file
    void Initialize ()
    {
        std::string SchemeString = "MUSCL";
        amrex::ParmParse pp("FiniteVolume");
        pp.query("Scheme",SchemeString);
        pp.query("ID_Bound",ID_Bound);
        pp.query("Quad",ID_Quad);
        if (SchemeString == "GODUNOV"){
            Scheme = FVReconstruction::GODUNOV;
            Num_Grow   = 1;
            ID_WENO    = 0;
            ID_Limiter = 0;
            ID_Poly    = 0;
        }
        else if (SchemeString == "MUSCL") {
            Scheme = FVReconstruction::MUSCL_MM;
            Num_Grow   = 2;
            ID_WENO    = 0;
            ID_Limiter = 1;
            ID_Poly    = 1;
            if (ID_Bound == 1){
                ID_Poly = 2;
            }
        }
        else if (SchemeString == "MUSCL-Minmod") {
            Scheme = FVReconstruction::MUSCL_MM;
            Num_Grow   = 2;
            ID_WENO    = 0;
            ID_Limiter = 1;
            ID_Poly    = 1;
            if (ID_Bound == 1){
                ID_Poly = 2;
            }
        }
        else if (SchemeString == "MUSCL-vanLeer") {
            Scheme = FVReconstruction::MUSCL_VL;
            Num_Grow   = 2;
            ID_WENO    = 0;
            ID_Limiter = 2;
            ID_Poly    = 1;
            if (ID_Bound == 1){
                ID_Poly = 2;
            }
        }
        else if (SchemeString == "MUSCL-Superbee") {
            Scheme = FVReconstruction::MUSCL_SB;
            Num_Grow   = 2;
            ID_WENO    = 0;
            ID_Limiter = 3;
            ID_Poly    = 1;
            if (ID_Bound == 1){
                ID_Poly = 2;
            }
        }
        else if (SchemeString == "MUSCL-MC") {
            Scheme = FVReconstruction::MUSCL_MC;
            Num_Grow   = 2;
            ID_WENO    = 0;
            ID_Limiter = 4;
            ID_Poly    = 1;
            if (ID_Bound == 1){
                ID_Poly = 2;
            }
        }
        else if (SchemeString == "WENO3") {
            Scheme = FVReconstruction::WENO_3;
            Num_Grow   = 2;
            ID_WENO    = 1;
            ID_Limiter = 3;
            ID_Poly    = 2;
        }
        else if (SchemeString == "WENO3B") {
            Scheme = FVReconstruction::WENO_3B;
            Num_Grow   = 2;
            ID_WENO    = 1;
            ID_Limiter = 4;
            ID_Poly    = 2;
        }
        else if (SchemeString == "WENO5") {
            Scheme = FVReconstruction::WENO_5;
            Num_Grow   = 3;
            ID_WENO    = 1;
            ID_Limiter = 5;
            ID_Poly    = 4;
        }
        else if (SchemeString == "WENO5B") {
            Scheme = FVReconstruction::WENO_5B;
            Num_Grow   = 3;
            ID_WENO    = 1;
            ID_Limiter = 6;
            ID_Poly    = 4;
        }
        else {
            amrex::Abort("FiniteVolume.Scheme = '" + SchemeString + "' is not recognized. Valid values are "
                         "GODUNOV, MUSCL, MUSCL-Minmod, MUSCL-vanLeer, MUSCL-Superbee, MUSCL-MC, "
                         "WENO3, WENO3B, WENO5, WENO5B.");
        }

        if (ID_Quad == 1){
#if (PHYSICS != FIVEEQS) || (AMREX_SPACEDIM != 2)
            amrex::Abort("FiniteVolume.Quad = 1 is only implemented for PHYSICS=FIVEEQS in 2D.");
#endif
            if (Scheme != FVReconstruction::WENO_5){
                amrex::Abort("FiniteVolume.Quad = 1 requires FiniteVolume.Scheme = WENO5.");
            }
            if (ID_Bound == 1){
                amrex::Abort("FiniteVolume.Quad = 1 cannot be combined with FiniteVolume.ID_Bound = 1.");
            }
        }

        //Quadrature_Weight
        if (ID_Poly == 0){
            Quadrature_Weight = 1.0;
        }
        else if (ID_Poly == 1){
            Quadrature_Weight = 0.5;
        }
        else if (ID_Poly == 2){
            Quadrature_Weight = 1.0/6.0;
        }
        else if (ID_Poly == 4){
            Quadrature_Weight = 1.0/12.0;
        }

        pp.query("THINC",ID_THINC);
        amrex::ParmParse pp_thinc("FiniteVolume.THINC");
        pp_thinc.query("UseNormal",ID_THINC_UseNormal);
        pp_thinc.query("EPS",THINC_EPS);
        pp_thinc.query("Beta",THINC_Beta);
        pp_thinc.query("NormalThreshold",THINC_NormalThreshold);
        pp_thinc.query("Multi",THINC_Multi);
        pp_thinc.query("FaceState",THINC_FaceState);
        pp_thinc.query("PhaseDensity",THINC_PhaseDensity);
        if (THINC_Multi < 0 || THINC_Multi > 1)               { amrex::Abort("FiniteVolume.THINC.Multi must be 0 or 1."); }
        if (THINC_FaceState < 0 || THINC_FaceState > 1)       { amrex::Abort("FiniteVolume.THINC.FaceState must be 0 or 1."); }
        if (THINC_PhaseDensity < 0 || THINC_PhaseDensity > 1) { amrex::Abort("FiniteVolume.THINC.PhaseDensity must be 0 or 1."); }
        ID_THINC_Options = (THINC_Multi != 0 || THINC_FaceState != 0 || THINC_PhaseDensity != 0) ? 1 : 0;
#if (PHYSICS != FIVEEQS) && (PHYSICS != FIVEEQS_NPHASE)
        if (ID_THINC_Options){
            amrex::Abort("FiniteVolume.THINC.Multi, FaceState and PhaseDensity are implemented for PHYSICS=FIVEEQS "
                         "and FIVEEQS_NPHASE only. Set them to 0 for this model.");
        }
#endif
        if (ID_THINC){
            // THINC reads two cells on each side of a face (GODUNOV has one ghost cell)
            Num_Grow = amrex::max(Num_Grow, 2);
            if (ID_Bound == 1){
                amrex::Print() << "Warning: FiniteVolume.THINC = 1 has no effect with FiniteVolume.ID_Bound = 1.\n";
            }
        }
    };
};

#endif