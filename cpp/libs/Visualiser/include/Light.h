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
#include "GLDebug.h"
#include "GLTypes.h"
#include "Shadow.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Light Abstract Base Class
***************************************************************************************************************************************************************/
enum class LightType
{
   Direct,
   Point,
   Spot,
   None
};

class Light
{
 protected:
   Light();

   Light(LightType type, const glm::vec4& rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity);

 public:
   Light(const Light& light) = delete;

   Light(Light&& light) noexcept = default;

   virtual ~Light() = default;

   virtual UInt Index() const = 0;

   virtual UInt LightCount() const = 0;

   inline void Init() { ShadowMap_.Init(2048, 2048); }

   inline auto& ShadowMap() { return ShadowMap_; }

   Light& operator=(const Light& light) = delete;

   Light& operator=(Light&& light) noexcept = default;

   inline void SetName(const std::string& name) { Name_ = name; }

   inline auto Name() const { return Name_; }

 protected:
   friend class Visualiser;
   friend class Shader;

   virtual void AddGUIElements();

   std::string Name_;
   glm::vec4   Colour_;
   GLfloat     AmbientIntensity_;
   GLfloat     DiffuseIntensity_;
   Shadow      ShadowMap_;
   LightType   Type_;
};

/***************************************************************************************************************************************************************
* Directional Light Class
***************************************************************************************************************************************************************/
class DirectLight : public Light
{
 public:
   DirectLight();

   DirectLight(glm::vec3 direction, glm::vec4 rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity);

   UInt Index() const override { return 0; }

   UInt LightCount() const override { return 1; }

   inline const auto& LightSpaceMatrix() const { return LightSpaceMatrix_; }

 private:
   friend class Shader;

   glm::vec3 Direction_;
   glm::mat4 LightSpaceMatrix_;
};

/***************************************************************************************************************************************************************
* Point Light Abstract Base Class
***************************************************************************************************************************************************************/
namespace detail {

template<class D>
class PointLightBase : public Light
{
 protected:
   PointLightBase() = delete;

   PointLightBase(LightType type, const glm::vec3& position, const glm::vec4& rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity,
                  const SVector3<GLfloat>& attenuation_coefficients);

 public:
   PointLightBase(const PointLightBase<D>& light) = delete;

   PointLightBase(PointLightBase<D>&& light) noexcept = default;

   ~PointLightBase();

   UInt Index() const override { return Index_; }

   UInt LightCount() const override { return PointLightCount_; }

   constexpr static auto FarPlane() { return FarPlane_; }

   inline const auto& Position() const { return Position_; }

   inline const auto& LightSpaceMatrices() const { return LightSpaceMatrices_; }

   PointLightBase<D>& operator=(const PointLightBase<D>& other) = delete;

   PointLightBase<D>& operator=(PointLightBase<D>&& other) = default;

 protected:
   void AddGUIElements() override;

   constexpr static UInt    MaxPointLights_{4};
   constexpr static GLfloat FarPlane_{25.0f};
   inline static UInt       PointLightCount_{0};
   UInt                     Index_;
   glm::vec3                Position_;
   SVector3<GLfloat>        AttenuationCoefficients_; // [0]: constant term, [1]: linear term, [2]: quadratic term
   SArray6<glm::mat4>       LightSpaceMatrices_;
};

}

/***************************************************************************************************************************************************************
* Point Light Class
***************************************************************************************************************************************************************/
class PointLight : public detail::PointLightBase<PointLight>
{
 public:
   PointLight(const glm::vec3& position, const glm::vec4& rgba_colour, GLfloat ambient_intensity, GLfloat diffuse_intensity,
              const SVector3<GLfloat>& attenuation_coefficients);

   void AddGUIElements() override;

   friend class Shader;
};

/***************************************************************************************************************************************************************
* Spotlight Class
***************************************************************************************************************************************************************/
class SpotLight : public detail::PointLightBase<SpotLight>
{
 public:
   SpotLight(const glm::vec3& position, const glm::vec3& direction, const glm::vec4& rgba_colour, GLfloat cone_angle, GLfloat ambient_intensity,
             GLfloat diffuse_intensity, const SVector3<GLfloat>& attenuation_coefficients);

   void AddGUIElements() override;

 private:
   friend class Shader;

   glm::vec3 Direction_;
   GLfloat   ConeAngle_;
   GLfloat   CosConeAngle_;
};

}
