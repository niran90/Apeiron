/***************************************************************************************************************************************************************
* GPL-3.0 License
* Copyright (C) 2022 Niran A. Ilangakoon
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program.
* If not, see <https://www.gnu.org/licenses/>.
***************************************************************************************************************************************************************/

#pragma once

namespace aprn {

namespace detail {

/** Macro expansions. */
#define EXPAND(x) x
#define PP_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define ARG_COUNT(...) EXPAND(PP_ARG_N(__VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))

/** Generic name definition for any macro with the number of arguments as its suffix. */
#define _MACRO_(name, n) name##n
#define _MACRO(name, n)  _MACRO_(name, n)
#define MACRO(func, ...) _MACRO(func, ARG_COUNT(__VA_ARGS__))(__VA_ARGS__)

}

/** Macros for defining constants. */
#define COMPILE_TIME_CONST(_constant, value) constexpr auto _constant = value;
#define RUN_TIME_CONST(_constant, value)     const auto _constant     = value;

}//aprn
