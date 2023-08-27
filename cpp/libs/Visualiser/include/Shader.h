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
#include "GLDebug.h"
#include "GLTypes.h"
#include "Camera.h"
#include "Light.h"
#include "Material.h"
#include "Texture.h"

#include <fstream>
#include <unordered_map>

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

enum class ShaderType
{
  Vertex,
  Geometry,
  Fragment,
  nTypes
};

struct ShaderSource
{
  std::string Vertex;
  std::string Geometry;
  std::string Fragment;
};

class Model;

class Shader
{
 public:
   Shader();

   Shader(const std::string& file_path);

   ~Shader();

   void Read(const std::string& file_path);

   inline void Bind() const { GLCall(glUseProgram(ID_)) }

   inline void Unbind() const { GLCall(glUseProgram(0)) }

   inline void SetWarnings(const bool is_on) { WarningsOn_ = is_on; }

   void UseModel(const Model& model);

   void UseMaterial(const Material& material);

   void UseTexture(const Texture& texture, const std::string& uniform_name, const UInt slot);

   void UseCamera(Camera& camera);

   void UseLight(const Light& light);

   void UsePostProcessor(bool status);

   void SetDirectionalShadowMap(const UInt slot);

   void SetPointShadowMap(const size_t i_point_light, const UInt slot);

   void SetPointPosition(const glm::vec3& position);

   void SetPointFarPlane(GLfloat far_plane);

   void SetDirectionalLightSpaceMatrix(const glm::mat4& light_space_matrix);

   void SetPointLightSpaceMatrices(const StaticArray<glm::mat4, 6>& light_space_matrices);

   void SetUniform1i(const std::string& name, GLint value);

   void SetUniform1f(const std::string& name, GLfloat value);

   void SetUniform2f(const std::string& name, GLfloat value0, GLfloat value1);

   void SetUniform3f(const std::string& name, GLfloat value0, GLfloat value1, GLfloat value2);

   void SetUniform4f(const std::string& name, GLfloat value0, GLfloat value1, GLfloat value2, GLfloat value3);

   void SetUniformMatrix4f(const std::string& name, const glm::mat4& proj_matrix);

   inline GLuint ID() const { return ID_; }

   inline static std::string Directory = "./libs/Visualiser/resources/shaders/";

private:
   ShaderSource Parse(const std::string& file_path);

   void Create(const std::string& vertex_shader, const std::string& geometry_shader, const std::string& fragment_shader);

   GLuint Compile(GLuint type, const std::string& source);

   void Attach(GLuint shader);

   void Delete(GLuint shader);

   void Delete();

   int UniformLocation(const std::string& name);

   GLuint                               ID_;
   std::unordered_map<std::string, int> UniformLocationCache_;
   bool                                 WarningsOn_{false};
};

}
