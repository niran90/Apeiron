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

#include "../include/Light.h"
#include "../include/GUI.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Light Abstract Base Class
***************************************************************************************************************************************************************/
Light::Light()
  : Light(LightType::None, {1.0f, 1.0f, 1.0f, 1.0f}, 1.0f, 0.0f) {}

Light::Light(LightType type, const glm::vec4& rgba_colour, const GLfloat ambient_intensity, const GLfloat diffuse_intensity)
  : Colour_(rgba_colour), AmbientIntensity_(ambient_intensity), DiffuseIntensity_(diffuse_intensity),
    ShadowMap_(type == OneOf(LightType::Point, LightType::Spot)), Type_(type) {}

void
Light::AddGUIElements()
{
   ImGui::ColorEdit4("Colour", &Colour_.r);
   ImGui::SliderFloat("Ambient Intensity", &AmbientIntensity_, 0.0f, 20.0f);
   ImGui::SliderFloat("Diffuse Intensity", &DiffuseIntensity_, 0.0f, 20.0f);
}

/***************************************************************************************************************************************************************
* Directional Light Class
***************************************************************************************************************************************************************/
DirectLight::DirectLight()
  : Light(), Direction_(glm::vec3(0.0, -1.0, 0.0)) {}

DirectLight::DirectLight(glm::vec3 direction, glm::vec4 rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity)
  : Light(LightType::Direct, rgba_colour, ambient_intensity, diffuse_intensity), Direction_(direction)
{
  const glm::mat4&& proj_matrix = glm::ortho(-5.0f, 5.0f, -5.0f, 5.0f, 1.0f, 20.0f);
  const glm::mat4&& view_matrix = glm::lookAt(-10.0f * Direction_, {0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
   LightSpaceMatrix_ = proj_matrix * view_matrix;
}

/***************************************************************************************************************************************************************
* Point Light Abstract Base Class
***************************************************************************************************************************************************************/
namespace detail {

template<class D>
PointLightBase<D>::PointLightBase(LightType type, const glm::vec3& position, const glm::vec4& rgba_colour,
                                        GLfloat ambient_intensity, GLfloat diffuse_intensity,
                                        const SVector3<GLfloat>& attenuation_coefficients)
  : Light(type, rgba_colour, ambient_intensity, diffuse_intensity), Position_(position), AttenuationCoefficients_(attenuation_coefficients)
{
  const auto proj_matrix = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, FarPlane_); // Note the aspect ratio of 1.0f

  // Create light space matrices for each cube map face.
  LightSpaceMatrices_[0] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3( 1.0, 0.0, 0.0), {0.0,-1.0, 0.0}); // Right face
  LightSpaceMatrices_[1] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3(-1.0, 0.0, 0.0), {0.0,-1.0, 0.0}); // Left face
  LightSpaceMatrices_[2] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3( 0.0, 1.0, 0.0), {0.0, 0.0, 1.0}); // Top face
  LightSpaceMatrices_[3] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3( 0.0,-1.0, 0.0), {0.0, 0.0,-1.0}); // Bottom face
  LightSpaceMatrices_[4] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3( 0.0, 0.0, 1.0), {0.0,-1.0, 0.0}); // Near face
  LightSpaceMatrices_[5] = proj_matrix * glm::lookAt(Position_, Position_ + glm::vec3( 0.0, 0.0,-1.0), {0.0,-1.0, 0.0}); // Far face

  Index_ = PointLightCount_++;
}

template<class D>
PointLightBase<D>::~PointLightBase() { PointLightCount_--; }

template<class D>
void
PointLightBase<D>::AddGUIElements()
{
   ImGui::SliderFloat3("Position"   , &Position_.x, -10.0f, 10.0f);
   ImGui::SliderFloat3("Attenuation", AttenuationCoefficients_.data(), 0.0f, 2.0f);
   Light::AddGUIElements();
}

}

template class detail::PointLightBase<PointLight>;
template class detail::PointLightBase<SpotLight>;

/***************************************************************************************************************************************************************
* Point Light Class
***************************************************************************************************************************************************************/
PointLight::PointLight(const glm::vec3& position, const glm::vec4& rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity,
                       const SVector3<GLfloat>& attenuation_coefficients)
  : PointLightBase(LightType::Point, position, rgba_colour, ambient_intensity, diffuse_intensity, attenuation_coefficients) {}

void
PointLight::AddGUIElements() { PointLightBase::AddGUIElements(); }

/***************************************************************************************************************************************************************
* Spotlight Class
***************************************************************************************************************************************************************/
SpotLight::SpotLight(const glm::vec3& position, const glm::vec3& direction, const glm::vec4& rgba_colour, GLfloat cone_angle, GLfloat ambient_intensity,
                     GLfloat diffuse_intensity, const SVector3<GLfloat>& attenuation_coefficients)
  : PointLightBase(LightType::Spot, position, rgba_colour, ambient_intensity, diffuse_intensity, attenuation_coefficients),
    Direction_(glm::normalize(direction)), ConeAngle_(cone_angle), CosConeAngle_(std::cos(DegToRad(cone_angle))) {}

void
SpotLight::AddGUIElements()
{
   PointLightBase::AddGUIElements();
   EXIT("Need to populate elements.")
}

}
