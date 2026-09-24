// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef LINALG_H_
#define LINALG_H_


namespace linalg
{



template <typename M>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
CrossProduct (M& C, M const & A, M const & B)
{
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
}

template <typename E>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Eye (E& EyeM, int NRow, int NCol)
{

    for (int j = 0; j < NCol; j++){
        for (int i = 0; i < NRow; i++){
            EyeM[i + j*NRow] = ((i == j) ? 1.0 : 0.0);
        }
    }
}


template <typename M, typename V, typename R>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void 
GEMV (R & b,
      M const & A,
      V const & x,
      int m,
      int n)
{

    for (int k = 0; k < m; k++){
        b[k] = 0.;
    }

    for (int k = 0; k < n; k++){
        for (int i = 0; i < m; i++){
            int l = i + m*k;
            b[i] += A[l]*x[k];
        }
    }


};


template <typename M,typename N,typename P>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void 
GEMM(M & C,
     N const & A,
     P const & B,
     int m,
     int n,
     int p)
{
    for (int k = 0; k < m*p; k++){
        C[k] = 0.;
    }

    for (int k = 0; k < n; k++){
        for (int i = 0; i < m; i++){
            for (int j = 0; j < p; j++){
                int la = i + m*k;
                int lb = k + n*j;
                int lc = i + m*j;
                C[lc] += A[la]*B[lb];
            }
        }
    }

};



template <typename T>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void 
transpose (T & AT, T & A, int n, int m)
{
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            int a = i + j*n;
            int at = j + i*m;
            AT[at] = A[a];
        }
    }
}



template <typename IM, typename M>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void 
invert(IM & inverse, M & matrix, int n)
{
    
}



}


#endif 

