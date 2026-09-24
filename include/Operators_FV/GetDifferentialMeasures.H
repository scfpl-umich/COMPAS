// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef GetDifferentialMeasures_H_
#define GetDifferentialMeasures_H_

template <typename TDX>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetVolume (amrex::Real & Vol,
           int ID_COORDSYS,
           AMREX_D_DECL(amrex::Real x10,
           amrex::Real x20,
           amrex::Real x30),
           TDX const & dx)
{
    if (ID_COORDSYS==0){ // Cartesian
        Vol = dx[0]*dx[1]*dx[2];
    }
    else if (ID_COORDSYS==1){ // Cylindrical
        Vol = 0.5*dx[0]*dx[1]*dx[2]*(2.*x10 + dx[0]);
    }
    else if (ID_COORDSYS==2){ // Spherical
#if (AMREX_SPACEDIM > 1)
        Vol = 0.333333333333333*(std::cos(x20) - std::cos(x20 + dx[1]))*dx[0]*dx[2]*(3.*x10*x10 + 3.*x10*dx[0] + dx[0]*dx[0]);
#endif
    }

}

template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetVolume_K (amrex::Real & Vol,
             int i,
             int j,
             int k,
             int ID_COORDSYS,
             TDX const & dx,
             TP const & prob_lo)
{
    AMREX_D_TERM(
    amrex::Real x10 = prob_lo[0] + i*dx[0];,
    amrex::Real x20 = prob_lo[1] + j*dx[1];,
    amrex::Real x30 = prob_lo[2] + k*dx[2];)

    GetVolume(Vol,ID_COORDSYS,AMREX_D_DECL(x10,x20,x30),dx);
}



template <typename TDX>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetTimeIntegral (amrex::Real & Integral,
                 int ID_COORDSYS,
                 AMREX_D_DECL(amrex::Real x10,
                 amrex::Real x20,
                 amrex::Real x30),
                 TDX const & dx)
{
    if (ID_COORDSYS==0){ // Cartesian
        Integral = dx[0]*dx[1]*dx[2];
    }
    else if (ID_COORDSYS==1){ // Cylindrical
        Integral = 0.333333333333333*dx[0]*dx[1]*dx[2]*(3.*x10*x10 + 3.*x10*dx[0] + dx[0]*dx[0]);
    }
    else if (ID_COORDSYS==2){ // Spherical
#if (AMREX_SPACEDIM > 1)
        Integral = 0.1*(std::pow(x10+dx[0],5.) - std::pow(x10,5.))*(std::cos(2.*x20 + dx[1])*std::sin(dx[1]) - dx[1])*dx[2];
#endif
    }

}


template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetTimeIntegral_K (amrex::Real & Integral,
                   int i,
                   int j,
                   int k,
                   int ID_COORDSYS,
                   TDX const & dx,
                   TP const & prob_lo)
{
    AMREX_D_TERM(
    amrex::Real x10 = prob_lo[0] + i*dx[0];,
    amrex::Real x20 = prob_lo[1] + j*dx[1];,
    amrex::Real x30 = prob_lo[2] + k*dx[2];)

    GetTimeIntegral(Integral,ID_COORDSYS,AMREX_D_DECL(x10,x20,x30),dx);
}



template <typename TDX>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetSurfaceIntegral (amrex::Real & Area,
                    int ID_COORDSYS,
                    int dir,
                    amrex::Real face,
                    AMREX_D_DECL(amrex::Real x10,
                    amrex::Real x20,
                    amrex::Real x30),
                    TDX const & dx)
{
    if (dir == 0){

        if (ID_COORDSYS==0){
            Area = dx[1]*dx[2];

        } else if (ID_COORDSYS==1){
            Area = (x10 + 0.5*(face + 1.)*dx[0])*dx[1]*dx[2];

        } if (ID_COORDSYS==2){
#if (AMREX_SPACEDIM > 1)
            Area = 0.25*(std::cos(x20) - std::cos(x20 + dx[1]))*std::pow(2.*x10 + (1. + face)*dx[0],2.)*dx[2];
#endif
        }

    }
#if (AMREX_SPACEDIM > 1)
    else if (dir == 1){

        if (ID_COORDSYS==0){
            Area = dx[0]*dx[2];

        } else if (ID_COORDSYS==1){
            Area = dx[0]*dx[2];
            
        } if (ID_COORDSYS==2){
            Area = 0.5*std::sin(x20 + 0.5*(1. + face)*dx[1])*dx[0]*dx[2]*(2.*x10 + dx[0]);
            
        }

    }
#endif
#if (AMREX_SPACEDIM > 2)
    else if (dir == 2){

        if (ID_COORDSYS==0){
            Area = dx[0]*dx[1];

        } else if (ID_COORDSYS==1){
            Area = 0.5*(2.*x10 + dx[0])*dx[0]*dx[1];
            
        } if (ID_COORDSYS==2){
            Area = 0.5*(2.*x10 + dx[0])*dx[0]*dx[1];
            
        }

    }
#endif
}


template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetSurfaceIntegral_K (amrex::Real & Area,
                      int i,
                      int j,
                      int k,
                      int ID_COORDSYS,
                      int dir,
                      amrex::Real face,
                      TDX const & dx,
                      TP const & prob_lo)
{
    AMREX_D_TERM(
    amrex::Real x10 = prob_lo[0] + i*dx[0];,
    amrex::Real x20 = prob_lo[1] + j*dx[1];,
    amrex::Real x30 = prob_lo[2] + k*dx[2];)

    GetSurfaceIntegral(Area,ID_COORDSYS,dir,face,AMREX_D_DECL(x10,x20,x30),dx);

}


template <typename TDX>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetFluxIntegral (amrex::Real & Integral,
                 int ID_COORDSYS,
                 int dir,
                 amrex::Real face,
                 AMREX_D_DECL(amrex::Real x10,
                 amrex::Real x20,
                 amrex::Real x30),
                 TDX const & dx)
{
    if (dir == 0){

        if (ID_COORDSYS==0){
            Integral = dx[1]*dx[2];

        } else if (ID_COORDSYS==1){
            Integral = std::pow(x10 + 0.5*(face + 1.)*dx[0],2.)*dx[1]*dx[2];

        } if (ID_COORDSYS==2){
#if (AMREX_SPACEDIM > 1)
            Integral = 0.03125*dx[2]*(dx[1] - std::cos(2.*x20 + dx[1])*std::sin(dx[1]))*std::pow(2.*x10 + (1.+face)*dx[0],4.);
#endif
        }

    }
#if (AMREX_SPACEDIM > 1)
    else if (dir == 1){

        if (ID_COORDSYS==0){
            Integral = dx[0]*dx[2];

        } else if (ID_COORDSYS==1){
            Integral = dx[0]*dx[2];
            
        } if (ID_COORDSYS==2){
            Integral = 0.333333333333333*std::pow(std::sin(x20 + 0.5*(1. + face)*dx[1]),2.)*dx[0]*dx[2]*(3.*x10*x10 + 3.*x10*dx[0] + dx[0]*dx[0]);
            
        }

    }
#endif
#if (AMREX_SPACEDIM > 2)
    else if (dir == 2){

        if (ID_COORDSYS==0){
            Integral = dx[0]*dx[1];

        } else if (ID_COORDSYS==1){
            Integral = 0.333333333333333*(3.*x10*x10 + 3.*x10*dx[0] + dx[0]*dx[0])*dx[0]*dx[1];
            
        } if (ID_COORDSYS==2){
            Integral = 0.333333333333333*(3.*x10*x10 + 3.*x10*dx[0] + dx[0]*dx[0])*dx[0]*dx[1];
            
        }

    }
#endif
}


template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetFluxIntegral_K (amrex::Real & Integral,
                   int i,
                   int j,
                   int k,
                   int ID_COORDSYS,
                   int dir,
                   amrex::Real face,
                   TDX const & dx,
                   TP const & prob_lo)
{
    AMREX_D_TERM(
    amrex::Real x10 = prob_lo[0] + i*dx[0];,
    amrex::Real x20 = prob_lo[1] + j*dx[1];,
    amrex::Real x30 = prob_lo[2] + k*dx[2];)

    GetFluxIntegral(Integral,ID_COORDSYS,dir,face,AMREX_D_DECL(x10,x20,x30),dx);

}





template <typename TDX>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetSourceIntegral (amrex::Real & Integral,
                   int ID_COORDSYS,
                   AMREX_D_DECL(amrex::Real x10,
                   amrex::Real x20,
                   amrex::Real x30),
                   TDX const & dx)
{
    if (ID_COORDSYS==0){ // Cartesian
        Integral = dx[0]*dx[1]*dx[2];
    }
    else if (ID_COORDSYS==1){ // Cylindrical
        Integral = 0.5*dx[0]*dx[1]*dx[2]*(2.*x10 + dx[0]);
    }
    else if (ID_COORDSYS==2){ // Spherical
#if (AMREX_SPACEDIM > 1)
        Integral = 0.25*(std::cos(x20) - std::cos(x20 + dx[1]))*dx[2]*(std::pow(x10+dx[0],4.) - std::pow(x10,4.));
#endif
    }

}


template <typename TDX, typename TP>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
GetSourceIntegral_K (amrex::Real & Integral,
                     int i,
                     int j,
                     int k,
                     int ID_COORDSYS,
                     TDX const & dx,
                     TP const & prob_lo)
{
    AMREX_D_TERM(
    amrex::Real x10 = prob_lo[0] + i*dx[0];,
    amrex::Real x20 = prob_lo[1] + j*dx[1];,
    amrex::Real x30 = prob_lo[2] + k*dx[2];)

    // amrex::Print() << "x10 " << x10 << "\n";

    GetSourceIntegral(Integral,ID_COORDSYS,AMREX_D_DECL(x10,x20,x30),dx);

}













#endif