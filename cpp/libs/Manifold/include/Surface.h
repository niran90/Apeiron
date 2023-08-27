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

#include "LinearAlgebra/include/Vector.h"

namespace aprn::mnfld {

class Surface
{
   using Vector    = SVectorR<3>;
   using Parameter = SVectorR<2>;

 public:
   constexpr Surface();

   virtual constexpr Vector Point(const Parameter& params) const = 0;

   virtual constexpr Vector Tangent(const Parameter& params) const = 0;

   virtual constexpr Vector Normal(const Parameter& params) const = 0;

   constexpr Vector Bitangent(const Vector& normal, const Vector& tangent) const;
};

/***************************************************************************************************************************************************************
* Linear/Piecewise Linear Surfaces
***************************************************************************************************************************************************************/
//class Plane : public Surface
//{
// public:
//   Plane(const SVectorR3& unit_normal, const SVectorR3& point = {Zero, Zero, Zero});
//};

}
