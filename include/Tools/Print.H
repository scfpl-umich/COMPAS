// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PPRINT_H_
#define PPRINT_H_

#include <cstdarg>
#include <string>
#include <utility>  // std::forward


// template<typename ...Args>
// constexpr void 
// PPrint (Args&&... args) noexcept
// {
//    ((amrex::Print() << std::forward<Args>(args) << " "), ...);
//    amrex::Print() << "\n";
// }


// AMREX_FORCE_INLINE
// void 
// PPrint (const char* fmt...) // C-style "const char* fmt, ..." is also valid
// {
//     va_list args;
//     va_start(args, fmt);
 
//     while (*fmt != '\0')
//     {
//         if (*fmt == 'd')
//         {
//             int i = va_arg(args, int);
//             amrex::Print() << i << ' ';
//         }
//         else if (*fmt == 'c')
//         {
//             // note automatic conversion to integral type
//             int c = va_arg(args, int);
//             amrex::Print() << static_cast<char>(c) << ' ';
//         }
//         else if (*fmt == 's')
//         {
//             // note automatic conversion to integral type
//             std::string s = va_arg(args, int);
//             amrex::Print() << s << ' ';
//         }
//         else if (*fmt == 'f')
//         {
//             double d = va_arg(args, double);
//             amrex::Print() << d << ' ';
//         }
//         ++fmt;
//     }
//     amrex::Print() << '\n';
 
//     va_end(args);
// }

#endif