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

#include <functional>
#include <memory>
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

class ModelFactory
{
   using Vector = SVector3<float>;
   using Point  = SVector3<float>;

   ModelFactory() = default;

 public:
   /** 1D parts
   ************************************************************************************************************************************************************/
   static Model Segment(const Point& p0, const Point& p1, float line_width = 0.1);

   static Model Segment(const Point& p0, const Point& p1, const Vector& miter0, const Vector& miter1, float line_width = 0.1);

   static Model SegmentChain(const DArray<Point>& points, float line_width = 0.1);

   static Model Arc(float radius, float angle, float line_width = Zero);

   static Model Arc(float radius, float start_angle, float end_angle, float line_width = Zero);

   /** 2D parts
   ************************************************************************************************************************************************************/
   static Model Triangle(float length, float border_width = Zero);

   static Model Triangle(float length, float height, float apex_ratio, float border_width = Zero);

   static Model Triangle(const SArray3<Point>& points, float border_width = Zero);

   static Model Square(float length, float border_width = Zero);

   static Model Rectangle(float length, float height, float border_width = Zero);

   static Model Quadrilateral(const SArray4<Point>& points, float border_width = Zero);

   static Model Polygon(const DArray<Point>& points, float border_width = Zero);

   static Model Sector(float radius, float angle, float border_width = Zero);

   static Model Sector(float radius, float start_angle, float end_angle, float border_width = Zero);

   static Model Circle(float radius, float border_width = Zero);

   static Model Circle(float radius, float start_angle, float border_width = Zero);

   static Model Ellipse(float radius_x, float radius_y, float border_width = Zero);

   static Model ScreenQuad();

   /** 3D parts
   ************************************************************************************************************************************************************/
   static Model Tetrahedron(float length);

   static Model Tetrahedron(const SArray4<Point>& points);

   static Model Cube(float length);

   static Model Cuboid(float length, float width, float height);

   static Model Octahedron(float length);

   static Model Dodecahedron(float length);

   static Model Icosahedron(float length);

   static Model Sphere(float radius);

   static Model Ellipsoid(float radius_x, float radius_y, float radius_z);

   static Model Cylinder(float radius, float height);

   static Model Cone(float radius, float height);
};

}
