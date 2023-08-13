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
#include "../../Manifold/include/Curve.h"
#include "GLDebug.h"
#include "GLTypes.h"
#include "Model.h"
#include "Object.h"

#include <functional>
#include <memory>
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

class ObjectFactory
{
   using Vector = SVector3<float>;
   using Point  = SVector3<float>;

   ObjectFactory() = default;

 public:
   /** 1D parts
   ************************************************************************************************************************************************************/
   static SPtr<Object> Segment(const Point& p0, const Point& p1, float line_width = 0.1);

   static SPtr<Object> Segment(const Point& p0, const Point& p1, const Vector& miter0, const Vector& miter1, float line_width = 0.1);

   static SPtr<Object> SegmentChain(const DArray<Point>& points, float line_width = 0.1);

   static SPtr<Object> Arc(float radius, float angle, float line_width = Zero);

   static SPtr<Object> Arc(float radius, float start_angle, float end_angle, float line_width = Zero);

   /** 2D parts
   ************************************************************************************************************************************************************/
   static SPtr<Object> Triangle(float length, float border_width = Zero);

   static SPtr<Object> Triangle(float length, float height, float apex_ratio, float border_width = Zero);

   static SPtr<Object> Triangle(const SArray3<Point>& points, float border_width = Zero);

   static SPtr<Object> Square(float length, float border_width = Zero);

   static SPtr<Object> Rectangle(float length, float height, float border_width = Zero);

   static SPtr<Object> Quadrilateral(const SArray4<Point>& points, float border_width = Zero);

   static SPtr<Object> Polygon(const DArray<Point>& points, float border_width = Zero);

   static SPtr<Object> Sector(float radius, float angle, float border_width = Zero);

   static SPtr<Object> Sector(float radius, float start_angle, float end_angle, float border_width = Zero);

   static SPtr<Object> Circle(float radius, float border_width = Zero);

   static SPtr<Object> Circle(float radius, float start_angle, float border_width = Zero);

   static SPtr<Object> Ellipse(float radius_x, float radius_y, float border_width = Zero);

   static SPtr<Object> ScreenQuad();

   /** 3D parts
   ************************************************************************************************************************************************************/
   static SPtr<Object> Tetrahedron(float length);

   static SPtr<Object> Tetrahedron(const SArray4<Point>& points);

   static SPtr<Object> Cube(float length);

   static SPtr<Object> Cuboid(float length, float width, float height);

   static SPtr<Object> Octahedron(float length);

   static SPtr<Object> Dodecahedron(float length);

   static SPtr<Object> Icosahedron(float length);

   static SPtr<Object> Sphere(float radius);

   static SPtr<Object> Ellipsoid(float radius_x, float radius_y, float radius_z);

   static SPtr<Object> Cylinder(float radius, float height);

   static SPtr<Object> Cone(float radius, float height);
};

}
