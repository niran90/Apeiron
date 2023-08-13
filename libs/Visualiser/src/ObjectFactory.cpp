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

#include "../include/ObjectFactory.h"

namespace aprn::vis {

/** 1D models
***************************************************************************************************************************************************************/
SPtr<Object>
ObjectFactory::Segment(const Point& p0, const Point& p1, const float line_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::Segment(const Point& p0, const Point& p1, const Vector& miter0, const Vector& miter1, const float line_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::SegmentChain(const DArray<Point>& points, const float line_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::Arc(const float radius, const float angle, const float line_width) { return Arc(radius, Zero, angle, line_width); }

SPtr<Object>
ObjectFactory::Arc(float radius, float start_angle, float end_angle, float line_width)
{
   return std::make_shared<Model>(Model());
}

/** 2D models
***************************************************************************************************************************************************************/
SPtr<Object>
ObjectFactory::Triangle(const float length, const float border_width) { return Triangle(length, length * std::sin(ThirdPi), Half, border_width); }

SPtr<Object>
ObjectFactory::Triangle(const float length, const float height, const float apex_ratio, const float border_width)
{
   const float x = Half * static_cast<float>(length);
   const float y = static_cast<float>(height);
   const float apex_x = x * (Two * apex_ratio - One);

   return Triangle({Point{-x, 0.0f, 0.0f}, Point{x, 0.0f, 0.0f}, Point{apex_x, y, 0.0f}}, border_width);
}

SPtr<Object>
ObjectFactory::Triangle(const SArray3<Point>& points, const float border_width) { return Polygon({points[0], points[1], points[2]}, border_width); }

SPtr<Object>
ObjectFactory::Square(const float length, const float border_width) { return Rectangle(length, length, border_width); }

SPtr<Object>
ObjectFactory::Rectangle(const float length, const float height, const float border_width)
{
   const float x = Half * length;
   const float y = Half * height;

   return Quadrilateral({Point{-x, -y, 0.0f}, Point{x, -y, 0.0f}, Point{x,  y, 0.0f}, Point{-x,  y, 0.0f}}, border_width);
}

SPtr<Object>
ObjectFactory::Quadrilateral(const SArray4<Point>& points, const float border_width)
{
   auto quad = Polygon(DArray<Point>(points.begin(), points.end()), border_width);
   auto& mesh = std::dynamic_pointer_cast<Model>(quad)->Mesh_;

   // Set tangents
   FOR_EACH(vertex, mesh.Vertices_) vertex.Tangent = glm::vec3(1.0f, 0.0f, 0.0f);

   // Set texture coordinates
   mesh.Vertices_[0].TextureCoordinates = glm::vec2(0.0f, 0.0f);
   mesh.Vertices_[1].TextureCoordinates = glm::vec2(1.0f, 0.0f);
   mesh.Vertices_[2].TextureCoordinates = glm::vec2(1.0f, 1.0f);
   mesh.Vertices_[3].TextureCoordinates = glm::vec2(0.0f, 1.0f);

   return quad;
}

SPtr<Object>
ObjectFactory::Polygon(const DArray<Point>& points, const float border_width)
{
   Model poly;
   poly.Mesh_.Shading_ = ShadingType::Flat;

   auto& vertices = poly.Mesh_.Vertices_;
   auto& indices  = poly.Mesh_.Indices_;

   const auto n_points = points.size();
   vertices.resize(n_points);
   FOR(i, n_points) vertices[i].Position = SVectorToGlmVec(points[i]);

   // TODO - extend to non-convex polygons as well!
   indices.resize(3 * (n_points - 2));
   FOR(i, n_points - 2)
   {
      indices[3 * i] = 0;
      FOR(j, 1, 3) indices[j + 3 * i] = i + j;
   }

   return std::make_shared<Model>(poly);
}

SPtr<Object>
ObjectFactory::Sector(const float radius, const float angle, const float border_width) { return Sector(radius, Zero, angle, border_width); }

SPtr<Object>
ObjectFactory::Sector(float radius, float start_angle, float end_angle, float border_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::Circle(const float radius, const float border_width) { return Circle(radius, Zero, border_width); }

SPtr<Object>
ObjectFactory::Circle(float radius, float start_angle, float border_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::Ellipse(const float radius_x, const float radius_y, const float border_width)
{
   return std::make_shared<Model>(Model());
}

SPtr<Object>
ObjectFactory::ScreenQuad() { return Square(Two, false); }

/** 3D models
***************************************************************************************************************************************************************/
SPtr<Object>
ObjectFactory::Tetrahedron(const float length)
{
   const float width  = static_cast<float>(length) * std::sin(ThirdPi);
   const float height = std::sqrt(iPow(length, 2) - iPow(TwoThird * width, 2));

   return Tetrahedron({Point{-0.5f * length, 0.0f, width / 3.0f},
                       Point{0.5f * length, 0.0f, width / 3.0f},
                       Point{0.0f, 0.0f, -2.0f * width / 3.0f},
                       Point{0.0f, height, 0.0f}});
}

SPtr<Object>
ObjectFactory::Tetrahedron(const SArray4<Point>& points)
{
   Model part;
   part.Mesh_.Shading_ = ShadingType::Flat;

   auto& vertices = part.Mesh_.Vertices_;
   auto& indices  = part.Mesh_.Indices_;

   vertices.resize(4);
   indices.resize(12);
   FOR(i, 4) vertices[i].Position = SVectorToGlmVec(points[i]);

   // Face 0
   indices[0] = 0;
   indices[1] = 2;
   indices[2] = 1;

   // Face 1
   indices[3] = 0;
   indices[4] = 1;
   indices[5] = 3;

   // Face 2
   indices[6] = 0;
   indices[7] = 3;
   indices[8] = 2;

   // Face 3
   indices[9]  = 1;
   indices[10] = 2;
   indices[11] = 3;

   return std::make_shared<Model>(part);
}

SPtr<Object>
ObjectFactory::Cube(const float length) { return Cuboid(length, length, length); }

SPtr<Object>
ObjectFactory::Cuboid(const float length, const float width, const float height)
{
   const float x = Half * static_cast<float>(length);
   const float y = Half * static_cast<float>(width);
   const float z = Half * static_cast<float>(height);

   SVector3<float> v0{-x, -y,  z};
   SVector3<float> v1{ x, -y,  z};
   SVector3<float> v2{-x, -y, -z};
   SVector3<float> v3{ x, -y, -z};
   SVector3<float> v4{-x,  y,  z};
   SVector3<float> v5{ x,  y,  z};
   SVector3<float> v6{-x,  y, -z};
   SVector3<float> v7{ x,  y, -z};

   Model part;
   part.Mesh_.Shading_ = ShadingType::Flat;

   auto& vertices = part.Mesh_.Vertices_;
   auto& indices  = part.Mesh_.Indices_;

   vertices.resize(8);
   vertices[0].Position = SVectorToGlmVec(v0);
   vertices[1].Position = SVectorToGlmVec(v1);
   vertices[2].Position = SVectorToGlmVec(v2);
   vertices[3].Position = SVectorToGlmVec(v3);
   vertices[4].Position = SVectorToGlmVec(v4);
   vertices[5].Position = SVectorToGlmVec(v5);
   vertices[6].Position = SVectorToGlmVec(v6);
   vertices[7].Position = SVectorToGlmVec(v7);

   indices.resize(36);
   // Face 0
   indices[0] = 0;
   indices[1] = 1;
   indices[2] = 4;

   indices[3] = 4;
   indices[4] = 1;
   indices[5] = 5;

   // Face 1
   indices[6] = 3;
   indices[7] = 2;
   indices[8] = 6;

   indices[9]  = 3;
   indices[10] = 6;
   indices[11] = 7;

   // Face 2
   indices[12] = 1;
   indices[13] = 0;
   indices[14] = 2;

   indices[15] = 2;
   indices[16] = 3;
   indices[17] = 1;

   // Face 3
   indices[18] = 4;
   indices[19] = 5;
   indices[20] = 7;

   indices[21] = 4;
   indices[22] = 7;
   indices[23] = 6;

   // Face 4
   indices[24] = 1;
   indices[25] = 3;
   indices[26] = 7;

   indices[27] = 1;
   indices[28] = 7;
   indices[29] = 5;

   // Face 5
   indices[30] = 0;
   indices[31] = 6;
   indices[32] = 2;

   indices[33] = 0;
   indices[34] = 4;
   indices[35] = 6;

   return std::make_shared<Model>(part);
}

}
