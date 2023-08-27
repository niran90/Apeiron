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

#include "../include/Polygon.h"
#include "../../Functional/include/Explicit.h"

namespace aprn {
namespace ptope {

/***************************************************************************************************************************************************************
* Static/Dynamic Polygons
***************************************************************************************************************************************************************/
template<PolytopeCategory cat, size_t dim>
constexpr Polygon<cat, dim>::Polygon(const Real _side_length)
{

}

template<PolytopeCategory cat, size_t dim>
template<class... static_vector>
constexpr Polygon<cat, dim>::Polygon(const static_vector&... vertices)
{

}

template<size_t dim>
Polygon<PolytopeCategory::Arbitrary2D, dim>::Polygon(size_t _n_vertices, const Real radius)
{

}

template<size_t dim>
template<class... static_vector>
Polygon<PolytopeCategory::Arbitrary2D, dim>::Polygon(const static_vector&... vertices)
{

}

/***************************************************************************************************************************************************************
* Triangles
***************************************************************************************************************************************************************/
template<size_t dim>
constexpr Triangle<dim>::Triangle(const SVectorR<dim>& v0, const SVectorR<dim>& v1, const SVectorR<dim>& v2)
{

}

template<size_t dim>
constexpr Triangle<dim>::Triangle(const Real length, const Real height, const Real _apex_ratio)
{

}

/***************************************************************************************************************************************************************
* Quadrilaterals
***************************************************************************************************************************************************************/
template<size_t dim>
constexpr Quadrilateral<dim>::Quadrilateral(const SVectorR<dim>& v0, const SVectorR<dim>& v1, const SVectorR<dim>& v2, const SVectorR<dim>& v3)
{

}

template<size_t dim>
constexpr Quadrilateral<dim>::Quadrilateral(const Real length, const Real height, const Real angle)
{

}

}
}
