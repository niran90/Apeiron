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
#include "../../DataContainer/include/Array.h"
#include "../../LinearAlgebra/include/Vector.h"
#include "Categories.h"

namespace aprn::ptope {

/***************************************************************************************************************************************************************
* Polytope Abstract Base Class
***************************************************************************************************************************************************************/
template<class D, PolytopeCategory cat, size_t dim>
struct Polytope
{
   constexpr virtual ~Polytope() = 0;

   /** Vertices/Faces Access */
   constexpr auto& Vertices() { return Derived().Vertices_; }

   constexpr const auto& Vertices() const { return Derived().Vertices_; }

   constexpr auto& Faces() { return Derived().Faces; }

   constexpr const auto& Faces() const { return Derived().Faces; }

 private:
   /** Derived Class Access */
   constexpr D& Derived() noexcept { return static_cast<D&>(*this); }

   constexpr const D& Derived() const noexcept { return static_cast<const D&>(*this); }
};

/***************************************************************************************************************************************************************
* Static Polytope Abstract Class
***************************************************************************************************************************************************************/
template<PolytopeCategory cat, size_t dim>
struct StaticPolytope : public Polytope<StaticPolytope<cat, dim>, cat, dim>
{
   constexpr StaticPolytope() { STATIC_ASSERT(isStaticPolytope<cat>(), "The polytope's information must be known at compile time.") }

   virtual constexpr ~StaticPolytope() = 0;

   StaticArray<SVectorR<dim>, PolytopeVertexCount<cat>()> Vertices;
   constexpr static auto Faces{GetPolytopeFaces<cat, dim>()};
};

/***************************************************************************************************************************************************************
* Dynamic Polytope Abstract Class
***************************************************************************************************************************************************************/
template<PolytopeCategory cat, size_t dim>
struct DynamicPolytope : public Polytope<DynamicPolytope<cat, dim>, cat, dim>
{
   DynamicPolytope() { STATIC_ASSERT(isDynamicPolytope<cat>(), "Rather use StaticPolytope if the polytope's information is known at compile-time.") }

   virtual ~DynamicPolytope() = 0;

   DynamicArray<SVectorR<dim>> Vertices;
   DynamicArray<DynamicArray<size_t>> Faces;
};

}