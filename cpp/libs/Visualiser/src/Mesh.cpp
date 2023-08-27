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

#include "../include/Mesh.h"
#include "../include/GLTypes.h"

namespace aprn::vis {

template<typename T>
void VertexAttributeLayout::AddAttribute(const GLuint n_values)
{
  constexpr GLenum type_enum = GLTypeEnum<T>();
  Attributes.emplace_back(type_enum, n_values, (type_enum == GL_UNSIGNED_BYTE ? GL_TRUE : GL_FALSE));
  Stride += n_values * GLTypeSize(type_enum);
}

Mesh::Mesh()
{
  VertexLayout_.AddAttribute<GLfloat>(3); // Position
  VertexLayout_.AddAttribute<GLfloat>(3); // Normal
  VertexLayout_.AddAttribute<GLfloat>(3); // Tangent
  VertexLayout_.AddAttribute<GLfloat>(4); // Colour
  VertexLayout_.AddAttribute<GLfloat>(2); // Texture coordinates

  DEBUG_ASSERT(VertexLayout_.Stride == sizeof(Vertex), "Check if the vertex attribute layout has changed.")
}

void
Mesh::ComputeVertexNormals()
{
  if(Shading_ == ShadingType::Flat)
  {
    for(size_t it = 0; it < Indices_.size(); it += 3)
    {
      // Assign the normal of the current triangle to each vertex normal.
      const GLuint iv0 = Indices_[it];
      const GLuint iv1 = Indices_[it + 1];
      const GLuint iv2 = Indices_[it + 2];
      const glm::vec3 face_normal = glm::cross(Vertices_[iv1].Position - Vertices_[iv0].Position, Vertices_[iv2].Position - Vertices_[iv0].Position);
      FOR(iv, 3) Vertices_[Indices_[it + iv]].Normal = face_normal;
    }
  }
  else if(Shading_ == ShadingType::Phong)
  {
    // Zero all vertex normals.
    FOR_EACH(vertex, Vertices_) vertex.Normal = glm::vec3(0.0, 0.0, 0.0);

    for(size_t it = 0; it < Indices_.size(); it += 3)
    {
      // Compute the normal of the current triangle and update the normal at each vertex in the current triangle
      const GLuint iv0 = Indices_[it];
      const GLuint iv1 = Indices_[it + 1];
      const GLuint iv2 = Indices_[it + 2];
      const glm::vec3 face_normal = glm::cross(Vertices_[iv1].Position - Vertices_[iv0].Position, Vertices_[iv2].Position - Vertices_[iv0].Position);
      FOR(iv, 3) Vertices_[Indices_[it + iv]].Normal += face_normal;
    }

    // Normalise all vertex normals.
    FOR_EACH(vertex, Vertices_) vertex.Normal = glm::normalize(vertex.Normal);
  }
  else EXIT("Unrecognised shading type prescribed.")

  if(Shading_ == ShadingType::Phong) FOR_EACH(vertex, Vertices_) vertex.Normal = glm::vec3(0.0, 0.0, 0.0);

  for(size_t it = 0; it < Indices_.size(); it += 3)
  {
    // Compute the normal of the current triangular face and update the normal at each vertex in the current triangle accordingly.
    const GLuint iv0 = Indices_[it];
    const GLuint iv1 = Indices_[it + 1];
    const GLuint iv2 = Indices_[it + 2];
    glm::vec3 face_normal = glm::cross(Vertices_[iv1].Position - Vertices_[iv0].Position, Vertices_[iv2].Position - Vertices_[iv0].Position);

    if(Shading_ == ShadingType::Flat)
    {
      face_normal = glm::normalize(face_normal);
      FOR(iv, 3) Vertices_[Indices_[it + iv]].Normal = face_normal;
    }
    else if(Shading_ == ShadingType::Phong) FOR(iv, 3) Vertices_[Indices_[it + iv]].Normal += face_normal;
    else EXIT("Unrecognised shading type prescribed.")
  }

  if(Shading_ == ShadingType::Phong) FOR_EACH(vertex, Vertices_) vertex.Normal = glm::normalize(vertex.Normal);
}

}

