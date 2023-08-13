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
#include "Animator.h"

#include <optional>
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

enum class ShadingType { Flat, Phong };

struct Vertex
{
   glm::vec3 Position;
   glm::vec3 Normal;
   glm::vec3 Tangent;
   glm::vec4 Colour{1.0f};
   glm::vec2 TextureCoordinates;
};

struct VertexAttribute
{
   VertexAttribute(const GLenum type, const GLuint n_components, const GLboolean is_normalised)
      : GLType(type), nComponents(n_components), isNormalised(is_normalised) {}

   GLenum    GLType;
   GLuint    nComponents;
   GLboolean isNormalised;
};

struct VertexAttributeLayout
{
   VertexAttributeLayout() : Stride(0) {};

   template<typename T>
   void AddAttribute(GLuint n_values);

   DArray<VertexAttribute> Attributes;
   GLuint                  Stride;
};

class Mesh
{
 public:
   Mesh();

   void ComputeVertexNormals();

   inline bool Loaded() const { return !Vertices_.empty(); }

   inline const auto& VertexLayout() const { return VertexLayout_; }

 private:
   friend class Model;
   friend class ObjectFactory;
   friend class TeXGlyph;

   VertexAttributeLayout VertexLayout_;
   DArray<Vertex>        Vertices_;
   DArray<GLuint>        Indices_;
   ShadingType           Shading_{ShadingType::Flat};
};

}
