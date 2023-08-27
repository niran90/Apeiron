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

namespace aprn {

/***************************************************************************************************************************************************************
* Cross/Dot/Outer Products
***************************************************************************************************************************************************************/
template<typename T, class D>
constexpr T
InnerProduct(const Vector<T, D>& vector0, const Vector<T, D>& vector1)
{
   const auto& v0 = vector0.Derived();
   const auto& v1 = vector1.Derived();
   return std::inner_product(v0.begin(), v0.end(), v1.begin(), static_cast<T>(Zero));
}

template<typename T, class D>
constexpr SVectorR3
CrossProduct(const Vector<T, D>& vector0, const Vector<T, D>& vector1)
{
   const auto& v0 = vector0.Derived();
   const auto& v1 = vector1.Derived();
   ASSERT(v0.size() == v1.size() && (v0.size() == 2 || v0.size() == 3), "Cross products can only be computed for 2D or 3D vectors.")

   return v0.size() == 2 ? SVectorR3{Zero, Zero, v0[0] * v1[1] - v0[1] * v1[0]} :
          SVectorR3{v0[1] * v1[2] - v0[2] * v1[1], v0[2] * v1[0] - v0[0] * v1[2], v0[0] * v1[1] - v0[1] * v1[0]};
}

/***************************************************************************************************************************************************************
* Lp-Norms
***************************************************************************************************************************************************************/
template<size_t p, typename T, class D>
requires Arithmetic<T>
constexpr Real
LpNorm(const Vector<T, D>& v)
{
   STATIC_ASSERT(0 < p && p < 4, "Lp-norms can currently only be computed for p = 1, 2, 3.")

   const auto& vector = v.Derived();
   switch(p)
   {
      case 1: return Sum(vector.begin(), vector.end());
      case 2: return std::sqrt(InnerProduct(vector, vector));
      case 3: throw("TODO");
   }
}

template<typename T, class D>
constexpr T
L1Norm(const Vector<T, D>& v) { return static_cast<T>(LpNorm<1>(v)); }

template<typename T, class D>
constexpr Real
L2Norm(const Vector<T, D>& v) { return LpNorm<2>(v); }

template<typename T, class D>
constexpr T
LInfNorm(const Vector<T, D>& v) { return MaxEntry(v.Derived().begin(), v.Derived().end()); }

template<typename T, class D>
constexpr Real
Magnitude(const Vector<T, D>& v) { return L2Norm(v); }

template<typename T, class D>
constexpr bool
isNormalised(const Vector<T, D>& v) { return isEqual(Magnitude(v), One); }

template<typename T, class D>
constexpr D
Normalise(const Vector<T, D>& v)
{
   const auto magn = Magnitude(v);
   return !isEqual(magn, Zero) ? v / magn : throw std::invalid_argument("Cannot normalise a vector of zero magnitude.");
}

/***************************************************************************************************************************************************************
* Vector Angle/Alignment/Rotation
***************************************************************************************************************************************************************/
template<bool orientangle = false, typename T, class D>
constexpr Real
ComputeAngle(const Vector<T, D>& v0, const Vector<T, D>& v1, const SVectorR3& orient = zAxis3)
{
   const auto smallangle = std::acos(InnerProduct(Normalise(v0), Normalise(v1)));
   if constexpr(!orientangle) return smallangle;
   else return Sgn(InnerProduct(CrossProduct(v0, v1), orient)) * smallangle;
}

template<typename T, class D>
constexpr bool
isAligned(const Vector<T, D>& v0, const Vector<T, D>& v1, const Real angle_thresh = TwelfthPi)
{
   const auto angle = ComputeAngle(v0, v1);
   return isBounded(angle_thresh, Zero, HalfPi) ? angle < angle_thresh || angle > (Pi - angle_thresh) :
                                                 throw std::domain_error("Angle threshold is out of bounds.");
}

template<typename T, class D>
constexpr D
RotateAbout(const Vector<T, D>& vector, const Real& angle, const SVectorR3& axis = zAxis3)
{
   EXIT("TODO")
   return D{};
}

template<typename T, class D>
constexpr D
RotateTowards(const Vector<T, D>& vector, const Real& angle, const Vector<T, D>& reference)
{
   EXIT("TODO")
   return D{};
}

}
