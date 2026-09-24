// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Sign_H_
#define Sign_H_

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
T Sign(T const& x)
{
    if (x > T(0.0)){
        return T(1.0);
    }
    else if (x < T(0.0)){
        return T(-1.0);
    }
    else{
        return T(0.0);
    }
}

template <typename T>
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
T SignNZ(T const& x)
{
    if (x > T(0.0)){
        return T(1.0);
    }
    else {
        return T(-1.0);
    }
}

#endif