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
#include "Categories.h"
#include "Polytope.h"

namespace aprn::ptope {

/***************************************************************************************************************************************************************
* Static/Dynamic Polyhedra
***************************************************************************************************************************************************************/
template<PolytopeCategory cat>
struct Polyhedron : public StaticPolytope<cat, 3>
{
  Polyhedron() { STATIC_ASSERT((isNPolytope<cat, 3>()), "A polyhedron must be a 3-polytope.") }

  Polyhedron(const Real _side_length);

  template<class... t_static_vector>
  Polyhedron(const t_static_vector&... vertices);
};

template<>
struct Polyhedron<PolytopeCategory::Arbitrary3D> : public DynamicPolytope<PolytopeCategory::Arbitrary3D, 3>
{
  Polyhedron();
};

/***************************************************************************************************************************************************************
* Tetrahedra
***************************************************************************************************************************************************************/
struct Tetrahedron : public Polyhedron<PolytopeCategory::Tetrahedron>
{
  Tetrahedron(const SVectorR3& v0, const SVectorR3& v1, const SVectorR3& v2, const SVectorR3& v3);
};

struct RegularTetrahedron : public Tetrahedron
{
  RegularTetrahedron(const Real _side_length);
};

struct TrirectangularTetrahedron : public Tetrahedron
{
  TrirectangularTetrahedron(const Real length, const Real height, const Real width, const bool _flip_horizontally = false);
};

/***************************************************************************************************************************************************************
* Static/Dynamic Polygons
***************************************************************************************************************************************************************/
//void CreateCube(Model &model, GLfloat _side_length);
//
//void CreateCuboid(Model &model, const SVectorF3& _lengths);
//
//void CreateOctahedron(Model &model, GLfloat _side_length);
//
//void CreateDodecahedron(Model &model, GLfloat _side_length);
//
//void CreateIcosahedron(Model &model, GLfloat _side_length);
//
//void CreateSphere(Model &model, GLfloat radius);
//
//void CreateEllipsoid(Model &model, GLfloat radius_x, GLfloat radius_y, GLfloat radius_z);
//
//void CreateCylinder(Model &model);
//
//void CreateCone(Model &model);

}
