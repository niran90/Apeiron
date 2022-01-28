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

namespace Apeiron {

/***************************************************************************************************************************************************************
* Light Abstract Base Class
***************************************************************************************************************************************************************/
enum class LightType
{
  Directional,
  Point,
  Spot,
  None
};

class Light
{
  friend class Shader;

protected:
  Light();

  Light(LightType _light_type, glm::vec4 _rgba_colour, GLfloat _ambient_intensity, GLfloat _diffuse_intensity);

public:
  ~Light() = default;

  virtual UInt GetIndex() const = 0;

  virtual UInt GetLightCount() const = 0;

  const Shadow& GetShadowMap() { return ShadowMap; }

protected:
  LightType Type;
  glm::vec4 Colour;
  GLfloat AmbientIntensity;
  GLfloat DiffuseIntensity;

  Shadow ShadowMap;
};

/***************************************************************************************************************************************************************
* Directional Light Class
***************************************************************************************************************************************************************/
class DirectionalLight : public Light
{
  friend class Shader;

public:
  DirectionalLight();

  DirectionalLight(glm::vec3 _direction, glm::vec4 _rgba_colour, GLfloat _ambient_intensity, GLfloat _diffuse_intensity);

  ~DirectionalLight() = default;

  UInt GetIndex() const override { return 0; }

  UInt GetLightCount() const override { return 1; }

  const glm::mat4& GetLightSpaceMatrix() const { return LightSpaceMatrix; }

private:
  glm::vec3 Direction;
  glm::mat4 LightSpaceMatrix;
};

/***************************************************************************************************************************************************************
* Point Light Abstract Base Class
***************************************************************************************************************************************************************/
namespace Detail {

template<class derived>
class PointLightBase : public Light
{
protected:
  PointLightBase() = delete;

  PointLightBase(LightType _light_type, const glm::vec3& _position, const glm::vec4& _rgba_colour, GLfloat _ambient_intensity, GLfloat _diffuse_intensity,
                 const SVector3<GLfloat>& _attenuation_coefficients);

  PointLightBase(const PointLightBase<derived>& _point_light_base);

public:
  ~PointLightBase();

  UInt GetIndex() const override { return iPointLight; }

  UInt GetLightCount() const override { return nPointLights; }

  constexpr static GLfloat GetFarPlane() { return FarPlane; }

  const glm::vec3& GetPosition() { return Position; }

  const StaticArray<glm::mat4, 6>& GetLightSpaceMatrices() const { return LightSpaceMatrices; }

  PointLightBase<derived>& operator=(const PointLightBase<derived>& _other) = delete;

  const PointLightBase<derived>& operator=(const PointLightBase<derived>& _other) const = delete;

protected:
  constexpr static UInt MaxPointLights{4};
  constexpr static GLfloat FarPlane{25.0f};

  static UInt nPointLights;
  UInt iPointLight;

  glm::vec3 Position;
  StaticVector<GLfloat, 3> AttenuationCoefficients; // 0: constant term, 1: linear term, 2: quadratic term
  StaticArray<glm::mat4, 6> LightSpaceMatrices;
};

}

/***************************************************************************************************************************************************************
* Point Light Class
***************************************************************************************************************************************************************/
class PointLight : public Detail::PointLightBase<PointLight>
{
  friend class Shader;

public:
  PointLight(const glm::vec3& _position, const glm::vec4& _rgba_colour, GLfloat _ambient_intensity, GLfloat _diffuse_intensity,
             const SVector3<GLfloat>& _attenuation_coefficients);
};

/***************************************************************************************************************************************************************
* Spotlight Class
***************************************************************************************************************************************************************/
class SpotLight : public Detail::PointLightBase<SpotLight>
{
  friend class Shader;

public:
  SpotLight(const glm::vec3& _position, const glm::vec3& _direction, const glm::vec4& _rgba_colour, GLfloat _cone_angle, GLfloat _ambient_intensity,
            GLfloat _diffuse_intensity, const SVector3<GLfloat>& _attenuation_coefficients);

private:
  glm::vec3 Direction;
  GLfloat ConeAngle;
  GLfloat CosConeAngle;
};

}
