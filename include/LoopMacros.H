// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef LOOPMACROS_H
#define LOOPMACROS_H

#include<AMReX.H>
#include<Macros.H>


#define D_LOOP_1(X) X
#define D_LOOP_2(X) D_LOOP_1(X), D_LOOP_1(X)
#define D_LOOP_3(X) D_LOOP_2(X), D_LOOP_1(X)
#define D_LOOP_4(X) D_LOOP_3(X), D_LOOP_1(X)
#define D_LOOP_5(X) D_LOOP_4(X), D_LOOP_1(X)
#define D_LOOP_6(X) D_LOOP_5(X), D_LOOP_1(X)
#define D_LOOP_7(X) D_LOOP_6(X), D_LOOP_1(X)
#define D_LOOP_8(X) D_LOOP_7(X), D_LOOP_1(X)
#define D_LOOP_9(X) D_LOOP_8(X), D_LOOP_1(X)
#define D_LOOP_10(X) D_LOOP_9(X), D_LOOP_1(X)
#define D_LOOP_11(X) D_LOOP_10(X), D_LOOP_1(X)
#define D_LOOP_12(X) D_LOOP_11(X), D_LOOP_1(X)
#define D_LOOP_13(X) D_LOOP_12(X), D_LOOP_1(X)
#define D_LOOP_14(X) D_LOOP_13(X), D_LOOP_1(X)
#define D_LOOP_15(X) D_LOOP_14(X), D_LOOP_1(X)
#define D_LOOP_16(X) D_LOOP_15(X), D_LOOP_1(X)
#define D_LOOP_17(X) D_LOOP_16(X), D_LOOP_1(X)
#define D_LOOP_18(X) D_LOOP_17(X), D_LOOP_1(X)
#define D_LOOP_19(X) D_LOOP_18(X), D_LOOP_1(X)
#define D_LOOP_20(X) D_LOOP_19(X), D_LOOP_1(X)
#define D_LOOP_21(X) D_LOOP_20(X), D_LOOP_1(X)
#define D_LOOP_22(X) D_LOOP_21(X), D_LOOP_1(X)
#define D_LOOP_23(X) D_LOOP_22(X), D_LOOP_1(X)
#define D_LOOP_24(X) D_LOOP_23(X), D_LOOP_1(X)

#define D_EXPAND(...) __VA_ARGS__

#define D_LOOP__(N, X) D_EXPAND(D_LOOP_ ## N)(X)
#define D_LOOP_(N, X) D_LOOP__(N, X)
#define D_LOOP(N, X) D_LOOP_(D_EXPAND(N), X)


#define T_LOOP_1(X) X
#define T_LOOP_2(X) T_LOOP_1(X) T_LOOP_1(X)
#define T_LOOP_3(X) T_LOOP_2(X) T_LOOP_1(X)
#define T_LOOP_4(X) T_LOOP_3(X) T_LOOP_1(X)
#define T_LOOP_5(X) T_LOOP_4(X) T_LOOP_1(X)
#define T_LOOP_6(X) T_LOOP_5(X) T_LOOP_1(X)
#define T_LOOP_7(X) T_LOOP_6(X) T_LOOP_1(X)
#define T_LOOP_8(X) T_LOOP_7(X) T_LOOP_1(X)
#define T_LOOP_9(X) T_LOOP_8(X) T_LOOP_1(X)
#define T_LOOP_10(X) T_LOOP_9(X) T_LOOP_1(X)
#define T_LOOP_11(X) T_LOOP_10(X) T_LOOP_1(X)
#define T_LOOP_12(X) T_LOOP_11(X) T_LOOP_1(X)
#define T_LOOP_13(X) T_LOOP_12(X) T_LOOP_1(X)
#define T_LOOP_14(X) T_LOOP_13(X) T_LOOP_1(X)
#define T_LOOP_15(X) T_LOOP_14(X) T_LOOP_1(X)
#define T_LOOP_16(X) T_LOOP_15(X) T_LOOP_1(X)
#define T_LOOP_17(X) T_LOOP_16(X) T_LOOP_1(X)
#define T_LOOP_18(X) T_LOOP_17(X) T_LOOP_1(X)
#define T_LOOP_19(X) T_LOOP_18(X) T_LOOP_1(X)
#define T_LOOP_20(X) T_LOOP_19(X) T_LOOP_1(X)
#define T_LOOP_21(X) T_LOOP_20(X) T_LOOP_1(X)
#define T_LOOP_22(X) T_LOOP_21(X) T_LOOP_1(X)
#define T_LOOP_23(X) T_LOOP_22(X) T_LOOP_1(X)
#define T_LOOP_24(X) T_LOOP_23(X) T_LOOP_1(X)

#define T_EXPAND(...) __VA_ARGS__

#define T_LOOP__(N, X) T_EXPAND(T_LOOP_ ## N)(X)
#define T_LOOP_(N, X) T_LOOP__(N, X)
#define T_LOOP(N, X) T_LOOP_(T_EXPAND(N), X)





#endif