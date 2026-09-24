// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PhaseField_Parameter_H_
#define PhaseField_Parameter_H_

#include<AMReX_ParmParse.H>


enum class PhaseFieldMechanism : int {
    CAC,
    CAC_Adv,
    CAC_VanDerWaals,
    CAC_Adv_NormalFace,
    CDI,
    CDI_CompressionInterpolation,
    CDI_NormalFace,
    ACDI,
    ACDI_NormalFace,
    CH,
    CH_VanDerWaals,
    CH_Adv,
    CH_Adv_NormalFace
};

struct PhaseField_Parameter
{
    int phase_field         = 0;
    int max_step            = 1000000;
    int Switch_Step         = 0;
    amrex::Real Switch_Time = 0.0;
    amrex::Real cfl         = 0.4;

    PhaseFieldMechanism Mechanism = PhaseFieldMechanism::CAC;
    int ID_FluxBasedMechanism = 0;
    int ID_Mechanism          = 1;

    amrex::Real Lambda     = 0.0; //User-input value; Lambda <= 0 will use a default value Lambda = 1.0
    amrex::Real Eta        = 0.0; //User-input value; Eta <=0 will use the finest-grid size
    amrex::Real Eta_Factor = 2.0;
    amrex::Real Eta_Multiplier = 1.0;
    amrex::Real Mobility   = 1.0;

    int ID_DegenerateMobility = 0; //0: Constant mobility
    int ID_Potential          = 0; //0: Double-Well potential
    int ID_Weight             = 0; //0: Quadratic weight

    int ID_ExplicitRC     = 0; //0: Implicit RC;        1: Explicit RC
    int ID_Flux           = 0; //0: UD flux;            1: HLL flux
    int ID_FluxPhase2Flux = 1; //0: Phasic formulation; 1: Mixture formulation
    int ID_Formulation    = 1; //0: Base formulation;   1: Efficient formulation

    //Read parameters from the inputs file
    void Initialize ()
    {
        std::string MechanismString = "CAC-Adv";
        amrex::ParmParse pp("PhaseField");
        pp.query("phase_field",phase_field);
        pp.query("max_step",max_step);
        pp.query("Switch_Step",Switch_Step);
        pp.query("Switch_Time",Switch_Time);
        pp.query("cfl",cfl);
        pp.query("Mechanism",MechanismString);
        pp.query("Lambda",Lambda);
        pp.query("Eta",Eta);
        pp.query("Eta_Factor",Eta_Factor);
        pp.query("Eta_Multiplier",Eta_Multiplier);
        //pp.query("Mobility",Mobility);
        pp.query("ID_DegenerateMobility",ID_DegenerateMobility);
        pp.query("ID_Potential",ID_Potential);
        pp.query("ID_Weight",ID_Weight);
        pp.query("ID_ExplicitRC",ID_ExplicitRC);
        pp.query("ID_Flux",ID_Flux);
        pp.query("ID_FluxPhase2Flux",ID_FluxPhase2Flux);
        pp.query("ID_Formulation",ID_Formulation);
        if (MechanismString == "CAC-VanDerWaals"){
            Mechanism             = PhaseFieldMechanism::CAC;
            ID_FluxBasedMechanism = 0;
            ID_Mechanism          = 0;
        }
        else if (MechanismString == "CAC-Adv"){
            Mechanism             = PhaseFieldMechanism::CAC;
            ID_FluxBasedMechanism = 0;
            ID_Mechanism          = 1;
        }
        else if (MechanismString == "CAC-Adv-NormalFace"){
            Mechanism             = PhaseFieldMechanism::CAC;
            ID_FluxBasedMechanism = 0;
            ID_Mechanism          = 2;
        }
        else if (MechanismString == "CDI"){
            Mechanism             = PhaseFieldMechanism::CDI;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 0;
        }
        else if (MechanismString == "CDI-CompressionInterpolation"){
            Mechanism             = PhaseFieldMechanism::CDI;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 1;
        }
        else if (MechanismString == "CDI-NormalFace"){
            Mechanism             = PhaseFieldMechanism::CDI;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 2;
        }
        else if (MechanismString == "ACDI"){
            Mechanism             = PhaseFieldMechanism::ACDI;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 3;
        }
        else if (MechanismString == "ACDI-NormalFace"){
            Mechanism             = PhaseFieldMechanism::ACDI;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 4;
        }
        else if (MechanismString == "CH-VanDerWaals"){
            Mechanism             = PhaseFieldMechanism::CH;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 5;
        }
        else if (MechanismString == "CH-Adv"){
            Mechanism             = PhaseFieldMechanism::CH;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 6;
        }
        else if (MechanismString == "CH-Adv-NormalFace"){
            Mechanism             = PhaseFieldMechanism::CH;
            ID_FluxBasedMechanism = 1;
            ID_Mechanism          = 7;
        }
        else {
            amrex::Abort("PhaseField.Mechanism = '" + MechanismString + "' is not recognized. Valid values are "
                         "CAC-Adv, CAC-Adv-NormalFace, CAC-VanDerWaals, CDI, CDI-NormalFace, "
                         "CDI-CompressionInterpolation, ACDI, ACDI-NormalFace, CH-Adv, CH-Adv-NormalFace, "
                         "CH-VanDerWaals.");
        }
        if (ID_FluxBasedMechanism == 0 && ID_ExplicitRC == 1){
            amrex::Abort("non flux-based mechanism can not use the explicit RC approach");
        }
    };
};

#endif