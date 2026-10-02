// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_FaceVelocityGradient_D_H_
#define FDM_FaceVelocityGradient_D_H_

// =======================================
// Velocity and velocity gradient on the lower face of cell (i,j,k) in direction ID_Dimension
// (face i-1/2, j-1/2 or k-1/2), from cell values, for the viscous fluxes of any model.
//
// The stencil is the one of the five-equation diffusion kernels
// (FDM_Conservative2FluxDiffusion_K_Compressible5Eq.H): 1 ghost cell with the edge neighbors,
// 8 cells in 2D and 16 in 3D. FaceStencil::Offset gives the cell of each stencil point
// (FaceStencil::OffsetOnFace the cell to read for a given face direction), and
// V[b][s] is velocity component b at stencil point s. The arithmetic is that of
// FDM_FluxDiffusion_D_Compressible5Eq.H:
//   face velocity       u_f = (u_m + u_c)/2
//   normal derivative   (u_c - u_m)/dx_d
//   tangential          from the edge averages of 4 cells, e.g. on an x face
//                       [u(i,j+1) + u(i-1,j+1) - u(i,j-1) - u(i-1,j-1)]/(4 dy)
// Grad_Vel[a][b] is d u_b / d x_a. Only the entries that the stress on this face needs are
// computed (the normal row and column and the diagonal); the others are 0.
// =======================================

namespace FaceStencil
{
enum Point : int {
    c = 0, xm,
#if (AMREX_SPACEDIM > 1)
    xp, ym, yp, xm_ym, xm_yp, xp_ym,
#endif
#if (AMREX_SPACEDIM > 2)
    zm, zp, xm_zm, xm_zp, xp_zm, ym_zm, ym_zp, yp_zm,
#endif
    N
};

// Index offset (di,dj,dk) of stencil point s from cell (i,j,k)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void Offset (int s, int& di, int& dj, int& dk)
{
    di = 0; dj = 0; dk = 0;
    switch (s) {
    case xm:    di = -1;                   break;
#if (AMREX_SPACEDIM > 1)
    case xp:    di =  1;                   break;
    case ym:             dj = -1;          break;
    case yp:             dj =  1;          break;
    case xm_ym: di = -1; dj = -1;          break;
    case xm_yp: di = -1; dj =  1;          break;
    case xp_ym: di =  1; dj = -1;          break;
#endif
#if (AMREX_SPACEDIM > 2)
    case zm:                      dk = -1; break;
    case zp:                      dk =  1; break;
    case xm_zm: di = -1;          dk = -1; break;
    case xm_zp: di = -1;          dk =  1; break;
    case xp_zm: di =  1;          dk = -1; break;
    case ym_zm:          dj = -1; dk = -1; break;
    case ym_zp:          dj = -1; dk =  1; break;
    case yp_zm:          dj =  1; dk = -1; break;
#endif
    default: break;
    }
}

// Offset of stencil point s for the face in direction ID_Dimension. The points one cell above
// the face cell (+1 along the face normal: xp, xp_ym, xp_zm on an x face) are not used on that
// face, and on the last face of a box they lie two cells outside it, beyond the single ghost
// cell of FiniteVolume.Scheme = GODUNOV; they get the offset of c, so that no stencil read
// leaves the first ghost layer. The values used on the face are unchanged.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void OffsetOnFace (int s, int ID_Dimension, int& di, int& dj, int& dk)
{
    Offset(s, di, dj, dk);
    int const dn = (ID_Dimension == 0) ? di : ((ID_Dimension == 1) ? dj : dk);
    if (dn > 0) { di = 0; dj = 0; dk = 0; }
}

// The neighbor across the face: xm, ym or zm
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int Minus (int ID_Dimension)
{
#if (AMREX_SPACEDIM > 2)
    if (ID_Dimension == 2) { return zm; }
#endif
#if (AMREX_SPACEDIM > 1)
    if (ID_Dimension == 1) { return ym; }
#endif
    amrex::ignore_unused(ID_Dimension);
    return xm;
}
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_FaceVelocityGradient_D(amrex::Array<amrex::Real,AMREX_SPACEDIM>& Vel,
                                amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM>& Grad_Vel,
                                amrex::Array<amrex::Array<amrex::Real,FaceStencil::N>,AMREX_SPACEDIM> const& V,
                   AMREX_D_DECL(amrex::Real const& dx,
                                amrex::Real const& dy,
                                amrex::Real const& dz),
                                int ID_Dimension)
{
    using namespace FaceStencil;
    for (int a = 0; a < AMREX_SPACEDIM; ++a){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Grad_Vel[a][b] = 0.0;
        }
    }

    //x
    if (ID_Dimension == 0){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Interpolation2Face_Linear_D(Vel[b], V[b][xm], V[b][c]);
            FDM_GradientX_XFaceFromCenter_D(Grad_Vel[0][b], V[b][c], V[b][xm], dx);
        }
#if (AMREX_SPACEDIM > 1)
        // du/dy and dv/dy
        for (int b = 0; b < 2; ++b){
            FDM_GradientY_XFaceFromCenter_D(Grad_Vel[1][b],
                                            V[b][c],
                                            V[b][xm], V[b][ym], V[b][yp],
                                            V[b][xm_ym], V[b][xm_yp],
                                            dy);
        }
#endif
#if (AMREX_SPACEDIM > 2)
        // du/dz and dw/dz
        for (int b = 0; b < 3; b += 2){
            FDM_GradientZ_XFaceFromCenter_D(Grad_Vel[2][b],
                                            V[b][c],
                                            V[b][xm], V[b][zm], V[b][zp],
                                            V[b][xm_zm], V[b][xm_zp],
                                            dz);
        }
#endif
    }
#if (AMREX_SPACEDIM > 1)
    //y
    else if (ID_Dimension == 1){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Interpolation2Face_Linear_D(Vel[b], V[b][ym], V[b][c]);
            FDM_GradientY_YFaceFromCenter_D(Grad_Vel[1][b], V[b][c], V[b][ym], dy);
        }
        // du/dx and dv/dx
        for (int b = 0; b < 2; ++b){
            FDM_GradientX_YFaceFromCenter_D(Grad_Vel[0][b],
                                            V[b][c],
                                            V[b][xm], V[b][xp], V[b][ym],
                                            V[b][xm_ym], V[b][xp_ym],
                                            dx);
        }
#if (AMREX_SPACEDIM > 2)
        // dv/dz and dw/dz
        for (int b = 1; b < 3; ++b){
            FDM_GradientZ_YFaceFromCenter_D(Grad_Vel[2][b],
                                            V[b][c],
                                            V[b][ym], V[b][zm], V[b][zp],
                                            V[b][ym_zm], V[b][ym_zp],
                                            dz);
        }
#endif
    }
#endif
#if (AMREX_SPACEDIM > 2)
    //z
    else if (ID_Dimension == 2){
        for (int b = 0; b < AMREX_SPACEDIM; ++b){
            Interpolation2Face_Linear_D(Vel[b], V[b][zm], V[b][c]);
            FDM_GradientZ_ZFaceFromCenter_D(Grad_Vel[2][b], V[b][c], V[b][zm], dz);
        }
        // du/dx and dw/dx
        for (int b = 0; b < 3; b += 2){
            FDM_GradientX_ZFaceFromCenter_D(Grad_Vel[0][b],
                                            V[b][c],
                                            V[b][xm], V[b][xp], V[b][zm],
                                            V[b][xm_zm], V[b][xp_zm],
                                            dx);
        }
        // dv/dy and dw/dy
        for (int b = 1; b < 3; ++b){
            FDM_GradientY_ZFaceFromCenter_D(Grad_Vel[1][b],
                                            V[b][c],
                                            V[b][ym], V[b][yp], V[b][zm],
                                            V[b][ym_zm], V[b][yp_zm],
                                            dy);
        }
    }
#endif

}

#endif
