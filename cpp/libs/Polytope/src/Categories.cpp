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

#include "../include/Categories.h"
#include "../../DataContainer/include/Array.h"

namespace aprn {
namespace ptope {

template<PolytopeCategory cat, size_t dim>
constexpr auto
GetPolytopeFaces(const StaticPolytope<cat, dim>&)
{
  STATIC_ASSERT(isStaticPolytope<cat>(), "The face vertices can only be determined for static polytopes.")

  constexpr auto n_faces = PolytopeFaceCount<cat>();
  typedef StaticArray<size_t, PolytopeFaceVertexCount<cat>()> face;
  typedef StaticArray<face, n_faces> faces;

  if(isNPolytope<cat, 2>()) return faces([](size_t i){ return face{i, (i + 1) % n_faces}; });
  else
  {
    switch(cat)
    {
      case PolytopeCategory::Tetrahedron:  return faces{{1, 0, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
      case PolytopeCategory::Cuboid:       return faces{{0, 1, 5, 4}, {1, 3, 7, 5}, {3, 2, 6, 7}, {2, 0, 4, 6}, {1, 0, 2, 3}, {4, 5, 7, 6}};
      case PolytopeCategory::Octahedron:   return faces{{0, 1, 4}, {1, 3, 4}, {3, 2, 4}, {2, 0, 4}, {1, 0, 5}, {3, 1, 5}, {2, 3, 5}, {0, 2, 5}};
      case PolytopeCategory::Dodecahedron: return faces{};
      case PolytopeCategory::Icosahedron:  return faces{};
      default: throw "Cannot determine the face vertices for the given polytope.";
    }
  }
}

}
}