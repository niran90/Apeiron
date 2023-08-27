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

#include "Constants.h"
#include "Debug.h"
#include "Types.h"
#include "BasicMath.h"

namespace aprn {

constexpr Real RelativeTolerance = Small;
constexpr Real ZeroTolerance     = Small;

/** Floating-point toleranced/exact equality check. */
template<bool exact_precision = false>
constexpr bool isEqual(const Real& a, const Real& b, const Real& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{ return exact_precision ? a == b : Abs(a - b) < (Min(Abs(a), Abs(b)) == Zero ? zero_tol : relative_tol * Max(Abs(a), Abs(b))); }

/** Floating-point toleranced/exact 'is lesser than' check. */
template<bool exact_precision = false>
constexpr bool isLess(const Real& a, const Real& b, const Real& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{ return a < b && !isEqual<exact_precision>(a, b, relative_tol, zero_tol); }

/** Floating-point toleranced/exact 'is lesser than or equal to' check. */
template<bool exact_precision = false>
constexpr bool isLessEqual(const Real& a, const Real& b, const Real& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{ return a <= b || (!exact_precision && isEqual(a, b, relative_tol, zero_tol)); }

/** Floating-point toleranced/exact 'is larger than' check. */
template<bool exact_precision = false>
constexpr bool isLarger(const Real& a, const Real& b, const Real& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{ return !isLessEqual<exact_precision>(a, b, relative_tol, zero_tol); }

/** Floating-point toleranced/exact 'is lesser than or equal to' check. */
template<bool exact_precision = false>
constexpr bool isLargerEqual(const Real& a, const Real& b, const Real& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{ return !isLess<exact_precision>(a, b, relative_tol, zero_tol); }

/** Floating-point toleranced/exact boundedness check. */
template<bool is_left_incl = true, bool is_right_incl = false, bool exact_precision = false, typename T>
requires Arithmetic<T>
constexpr bool isBounded(const T& a, const T& min, const T& max, const T& relative_tol = RelativeTolerance, const Real& zero_tol = ZeroTolerance)
{
   if constexpr(!isTypeSame<T, Real>()) return (is_left_incl ? min <= a : min < a) && (is_right_incl ? a <= max : a < max);
   else return (is_left_incl  ? isLessEqual<exact_precision>(min, a, relative_tol, zero_tol) : isLess<exact_precision>(min, a, relative_tol, zero_tol)) &&
               (is_right_incl ? isLessEqual<exact_precision>(a, max, relative_tol, zero_tol) : isLess<exact_precision>(a, max, relative_tol, zero_tol));
}

}
