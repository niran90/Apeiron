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

#include "../../../include/Global.h"
#include "../../LinearAlgebra/include/Vector.h"

namespace aprn::func {

/***************************************************************************************************************************************************************
* Functions from R -> R
***************************************************************************************************************************************************************/
constexpr Real
Linear(const Real x, const Real c0, const Real c1) { return c0 + c1 * x; }

constexpr Real
Quadratic(const Real x, const Real c0, const Real c1, const Real c2) { return c0 + c1 * x + c2 * x * x; }

constexpr Real
Cubic(const Real x, const Real c0, const Real c1, const Real c2, const Real c3) { return c0 + c1 * x + c2 * x * x + c3 * x * x * x; }

/***************************************************************************************************************************************************************
* Functions from R -> R^n
***************************************************************************************************************************************************************/
constexpr SVectorR2
Ellipse(const SVectorR2& radii, const Real theta) { return {radii[0] * std::cos(theta), radii[1] * std::sin(theta)}; }

constexpr SVectorR2
Circle(const Real radius, const Real theta) { return Ellipse({radius, radius}, theta); }

constexpr SVectorR3
Ellipsoid(const SVectorR3& _radii, const Real _theta, const Real _phi)
{
   return { _radii[0] * std::cos(_theta) * std::sin(_phi), _radii[1] * std::sin(_theta) * std::sin(_phi), _radii[2] * std::cos(_phi) };
}

constexpr SVectorR3
Sphere(const Real radius, const Real _theta, const Real _phi) { return Ellipsoid({radius, radius, radius}, _theta, _phi); }

//template<int p, int dimension = 1, typename T>
//constexpr T Polynomial(const StaticArray<T, dimension>& point, )
//{
//  STATIC_ASSERT(0 <= degree, "A polynomial cannot have a negative degree")
//  STATIC_ASSERT(0 < dimension, "The dimension of the polynomial must be at least 1")
//}
//
//template<int degree, int dimension = 1, typename T>
//constexpr T Polynomial(const StaticArray<T, dimension>& point, )
//{
//  STATIC_ASSERT(0 <= degree, "A polynomial cannot have a negative degree")
//  STATIC_ASSERT(0 < dimension, "The dimension of the polynomial must be at least 1")
//}

}