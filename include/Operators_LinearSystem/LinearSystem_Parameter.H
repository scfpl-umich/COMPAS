// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef LinearSystem_Parameter_H_
#define LinearSystem_Parameter_H_

#include<AMReX_ParmParse.H>
#include<AMReX_MLLinOp.H>
#include<AMReX_MLMG.H>

struct LinearSystem_Parameter
{
    bool composite_solve = true;

    int LinOp_maxorder = 3;

    bool agglomeration = true;
    bool consolidation = true;
    bool semicoarsening = false;

    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;

    amrex::Real tol_rel = 1.0e-10;
    amrex::Real tol_abs = 0.0;

    int max_iter     = 100;
    int max_fmg_iter = 0;
    int fixed_iter   = 100;

    int verbose        = 2;
    int bottom_verbose = 0;

    amrex::BottomSolver bottom_solver = MLMG::BottomSolver::bicgstab;
    int ID_BottomSolver = 0;

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("LinearSystem");
        pp.query("composite_solve",composite_solve);
        pp.query("LinOp_maxorder",LinOp_maxorder);
        pp.query("agglomeration",agglomeration);
        pp.query("consolidation",consolidation);
        pp.query("semicoarsening",semicoarsening);
        pp.query("max_coarsening_level",max_coarsening_level);
        pp.query("max_semicoarsening_level",max_semicoarsening_level);
        pp.query("tol_rel",tol_rel);
        pp.query("tol_abs",tol_abs);
        pp.query("max_iter",max_iter);
        pp.query("max_fmg_iter",max_fmg_iter);
        pp.query("fixed_iter",fixed_iter);
        pp.query("verbose",verbose);
        pp.query("bottom_verbose",bottom_verbose);
        pp.query("ID_BottomSolver",ID_BottomSolver);

        if (ID_BottomSolver == 0){
            bottom_solver = MLMG::BottomSolver::bicgstab; //The default.
        }
        else if (ID_BottomSolver == 1){
            bottom_solver = MLMG::BottomSolver::smoother; //Smoother such as Gauss-Seidel.
        }
        else if (ID_BottomSolver == 2){
            bottom_solver = MLMG::BottomSolver::bicgcg; //Start with bicgstab. Switch to cg if bicgstab fails. The matrix must be symmetric.
        }
        else if (ID_BottomSolver == 3){
            bottom_solver = MLMG::BottomSolver::cgbicg; //Start with cg. Switch to bicgstab if cg fails. The matrix must be symmetric.
        }
        else if (ID_BottomSolver == 4){
            bottom_solver = MLMG::BottomSolver::cg; //The conjugate gradient method. The matrix must be symmetric.
        }
    };
};

#endif