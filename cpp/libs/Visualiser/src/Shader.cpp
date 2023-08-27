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

#include "../include/Model.h"
#include "../include/Shader.h"

#include <glm/gtc/type_ptr.hpp>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Shader Public Interface
***************************************************************************************************************************************************************/
Shader::Shader()
  : ID_(0)
{
   GLCall(ID_ = glCreateProgram())
   ASSERT(ID_, "Could not create shader program.")
}

Shader::Shader(const std::string& file_path)
   : Shader() { Read(file_path); }

Shader::~Shader() { Delete(); }

void
Shader::Read(const std::string &file_path)
{
   auto source = Parse(file_path);
   Create(source.Vertex, source.Geometry, source.Fragment);
}

void
Shader::UseModel(const Model& model) { SetUniformMatrix4f("u_model_matrix", model.ModelMatrix()); }

void
Shader::UseMaterial(const Material& material)
{
   SetUniform1f("u_material.SpecularIntensity", material.SpecularIntensity_);
   SetUniform1f("u_material.Smoothness"       , material.Smoothness_);
}

void
Shader::UseTexture(const Texture& texture, const std::string& uniform_name, const UInt slot)
{
   SetUniform1i(uniform_name, slot);
   texture.Bind(slot);
}

void
Shader::UseCamera(Camera& camera)
{
   camera.UpdateViewMatrix();
   camera.UpdateProjMatrix();

   SetUniformMatrix4f("u_view_matrix"      , camera.ViewMatrix());
   SetUniformMatrix4f("u_projection_matrix", camera.ProjMatrix());

   const auto& camera_position = camera.Position();
   SetUniform3f("u_camera_position", camera_position.x, camera_position.y, camera_position.z);
}

void
Shader::UseLight(const Light& light)
{
   const LightType type = light.Type_;
   const UInt id        = light.Index();

   std::string uniform_name   = type == LightType::Direct ? "u_direct_light" :
                                type == LightType::Point  ? "u_point_lights[" + ToString(id) + "]" :
                                type == LightType::Spot   ? "u_spot_lights[" + ToString(id) + "]" :
                                throw std::invalid_argument("The lighting type was either not recognised or not specified.");

   std::string light_position = type == LightType::Point  ? "u_point_light_positions[" + ToString(id) + "]" :
                                type == LightType::Spot   ? "u_spot_light_positions[" + ToString(id) + "]" :
                                type == LightType::Direct ? "\0" :
                                throw std::invalid_argument("The lighting type was either not recognised or not specified.");

   std::string base_suffix    = type == LightType::Spot ? ".Point" : "\0";

   SetUniform4f(uniform_name + base_suffix + ".Base.Colour", light.Colour_.r, light.Colour_.g, light.Colour_.b, light.Colour_.a);
   SetUniform1f(uniform_name + base_suffix + ".Base.AmbientIntensity", light.AmbientIntensity_);
   SetUniform1f(uniform_name + base_suffix + ".Base.DiffuseIntensity", light.DiffuseIntensity_);

   if(type == LightType::Direct)
   {
      const glm::vec3& direction = static_cast<const DirectLight&>(light).Direction_;
      SetUniform3f(uniform_name + ".Direction", direction.x, direction.y, direction.z);
   }
   else if(type == LightType::Point)
   {
      const auto& point_light = static_cast<const PointLight&>(light);
      const auto& position    = point_light.Position_;
      const auto& att_coeffs  = point_light.AttenuationCoefficients_;

      SetUniform1i("u_point_light_count", light.LightCount());
      SetUniform3f(light_position, position.x, position.y, position.z);
      SetUniform3f(uniform_name + ".AttenuationCoefficients", att_coeffs[0], att_coeffs[1], att_coeffs[2]);
   }
   else if(type == LightType::Spot)
   {
      const auto& spot_light = static_cast<const SpotLight&>(light);
      const auto& position   = spot_light.Position_;
      const auto& att_coeffs = spot_light.AttenuationCoefficients_;
      const auto& direction  = spot_light.Direction_;

      SetUniform1i("u_spot_light_count", light.LightCount());
      SetUniform3f(light_position, position.x, position.y, position.z);
      SetUniform3f(uniform_name + ".Point.AttenuationCoefficients", att_coeffs[0], att_coeffs[1], att_coeffs[2]);
      SetUniform3f(uniform_name + ".Direction", direction.x, direction.y, direction.z);
      SetUniform1f(uniform_name + ".CosConeAngle", spot_light.CosConeAngle_);
   }
   else EXIT("Cannot yet handle the given light type.")
}

void
Shader::UsePostProcessor(bool status) { SetUniform1i("u_use_postprocessor", static_cast<GLint>(status)); }

void Shader::SetDirectionalShadowMap(const UInt slot) { SetUniform1i("u_direct_light.Shadow", slot); }

void Shader::SetPointShadowMap(const size_t i_point_light, const UInt slot) { SetUniform1i("u_point_lights[" + ToString(i_point_light) + "].Shadow", slot); }

void Shader::SetPointPosition(const glm::vec3& position) { SetUniform3f("u_plight_position", position.x, position.y, position.z); }

void Shader::SetPointFarPlane(const GLfloat far_plane) { SetUniform1f("u_point_light_far_plane", far_plane); }

void Shader::SetDirectionalLightSpaceMatrix(const glm::mat4& light_space_matrix) { SetUniformMatrix4f("u_dlight_space_matrix", light_space_matrix); }

void Shader::SetPointLightSpaceMatrices(const StaticArray<glm::mat4, 6>& light_space_matrices)
{
   FOR(i, 6) SetUniformMatrix4f("u_plight_space_matrices[" + ToString(i) + "]", light_space_matrices[i]);
}

/***************************************************************************************************************************************************************
* Setting Shader Uniforms
***************************************************************************************************************************************************************/
void Shader::SetUniform1i(const std::string& name, const GLint value) { GLCall(glUniform1i(UniformLocation(name), value)) }

void Shader::SetUniform1f(const std::string& name, const GLfloat value) { GLCall(glUniform1f(UniformLocation(name), value)) }

void Shader::SetUniform2f(const std::string& name, const GLfloat value0, const GLfloat value1) { GLCall(glUniform2f(UniformLocation(name), value0, value1)) }

void Shader::SetUniform3f(const std::string& name, const GLfloat value0, const GLfloat value1, const GLfloat value2)
{
   GLCall(glUniform3f(UniformLocation(name), value0, value1, value2))
}

void Shader::SetUniform4f(const std::string& name, const GLfloat value0, const GLfloat value1, const GLfloat value2, const GLfloat value3)
{
   GLCall(glUniform4f(UniformLocation(name), value0, value1, value2, value3))
}

void Shader::SetUniformMatrix4f(const std::string& name, const glm::mat4& proj_matrix)
{
   GLCall(glUniformMatrix4fv(UniformLocation(name), 1, GL_FALSE, glm::value_ptr(proj_matrix)))
}

/***************************************************************************************************************************************************************
* Shader Private Interface
***************************************************************************************************************************************************************/
ShaderSource
Shader::Parse(const std::string& file_path)
{

   std::string       line;
   std::stringstream string_stream[static_cast<UInt>(ShaderType::nTypes)];
   std::ifstream     stream(file_path);
   ShaderType        type{};

   while(getline(stream, line))
   {
      if(line.find("#shader") != std::string::npos)
      {
         if     (line.find("vertex")   != std::string::npos) type = ShaderType::Vertex;
         else if(line.find("geometry") != std::string::npos) type = ShaderType::Geometry;
         else if(line.find("fragment") != std::string::npos) type = ShaderType::Fragment;
      }
      else string_stream[static_cast<int>(type)] << line << '\n';
   }

   return { string_stream[0].str(), string_stream[1].str(), string_stream[2].str() };
}

void
Shader::Create(const std::string& vertex_shader, const std::string& geometry_shader, const std::string& fragment_shader)
{
   // Create and attach vertex shader.
   GLuint vs_id = Compile(GL_VERTEX_SHADER, vertex_shader);
   Attach(vs_id);

   // Create and attach geometry shader if it exists.
   GLuint gs_id(-1);
   if(!geometry_shader.empty())
   {
      gs_id = Compile(GL_GEOMETRY_SHADER, geometry_shader);
      Attach(gs_id);
   }

   // Create and attach fragment shader.
   GLuint fs_id = Compile(GL_FRAGMENT_SHADER, fragment_shader);
   Attach(fs_id);

   // Link shader program
   GLint result = 0;
   GLchar error_log[1024] = {0};
   GLCall(glLinkProgram(ID_))
   GLCall(glGetProgramiv(ID_, GL_LINK_STATUS, &result))
   if(!result)
   {
     GLCall(glGetProgramInfoLog(ID_, sizeof(error_log), nullptr, error_log))
     EXIT("Could not link shader program.")
   }

   // Validate shader program
   GLCall(glValidateProgram(ID_))
   GLCall(glGetProgramiv(ID_, GL_VALIDATE_STATUS, &result))
   if(!result)
   {
      GLCall(glGetProgramInfoLog(ID_, sizeof(error_log), nullptr, error_log))
      EXIT("Could not validate shader program.")
   }

   // Delete all shaders
   Delete(vs_id);
   if(!geometry_shader.empty()) Delete(gs_id);
   Delete(fs_id);
}

UInt
Shader::Compile(unsigned int type, const std::string& source)
{
   GLuint shader_ID = glCreateShader(type);
   const GLchar* src       = source.c_str();

   GLCall(glShaderSource(shader_ID, 1, &src, nullptr))
   GLCall(glCompileShader(shader_ID))

   // Check that shader compiled
   GLint result = 0;
   GLchar error_log[1024] = { 0 };
   GLCall(glGetShaderiv(shader_ID, GL_COMPILE_STATUS, &result))
   if(!result)
   {
      GLCall(glGetShaderInfoLog(shader_ID, sizeof(error_log), nullptr, error_log))
      EXIT("Failed to compile ", (type == GL_VERTEX_SHADER ? "vertex" : type == GL_GEOMETRY_SHADER ? "geometry" : "fragment"), " shader:\n ", error_log)
   }

   return shader_ID;
}

void Shader::Attach(const GLuint shader) { GLCall(glAttachShader(ID_, shader)) }

void Shader::Delete(const GLuint shader) { GLCall(glDeleteShader(shader)) }

void Shader::Delete() { GLCall(glDeleteProgram(ID_)) }

int Shader::UniformLocation(const std::string& name)
{
   if(UniformLocationCache_.contains(name)) return UniformLocationCache_[name];

   int location = glGetUniformLocation(ID_, name.c_str());
   if(WarningsOn_ && location < 0) WARN("Could not find the location for uniform ", name)
   UniformLocationCache_[name] = location;

   return location;
}

}

