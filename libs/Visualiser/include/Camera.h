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

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/vector_angle.hpp>

namespace aprn::vis {

class Camera
{
  constexpr static GLint KeyCount = 1024;

 public:
   Camera();

   Camera(const glm::vec3& position, GLfloat pitch, GLfloat yaw);

   void SetOrientation(const glm::vec3& position, GLfloat pitch, GLfloat yaw);

   void SetViewFrustum(const GLfloat& aspect_ratio, const GLfloat& field_of_view = -1.0, const GLfloat& near_plane = 0.0, const GLfloat& far_plane = 0.0);

   void KeyControl(const StaticArray<Bool, KeyCount>& keys, const GLfloat& delta_time);

   void CursorControl(const SVectorF2& cursor_displacement);

   void WheelControl(const SVectorF2& wheel_displacement);

   void UpdateViewMatrix();

   void UpdateProjectionMatrix();

   inline const glm::vec3& Position() const { return _Position; }

   inline const glm::mat4& ViewMatrix() const { return _ViewMatrix; }

   inline const glm::mat4& ProjectionMatrix() const { return _ProjMatrix; }

 private:
   void ClipPitch() { Clip(_Pitch, -89.0f, 89.0f); }

   glm::mat4 _ViewMatrix{1.0};
   glm::mat4 _ProjMatrix{1.0};
   glm::vec3 _Position;
   glm::vec3 _Front;
   glm::vec3 _Up;
   glm::vec3 _Right;
   GLfloat   _Pitch;
   GLfloat   _Yaw;
   GLfloat   _AspectRatio;
   GLfloat   _FieldOfView;
   GLfloat   _NearPlane;
   GLfloat   _FarPlane;

   constexpr static glm::vec3 _WorldUp    {0.0f, 1.0f, 0.0f};
   constexpr static GLfloat   _MotionSensitivity      {6.0f};
   constexpr static GLfloat   _RotationSensitivity    {7.0f};
   constexpr static GLfloat   _MouseCursorSensitivity {0.4f};
   constexpr static GLfloat   _MouseWheelSensitivity  {1.0f};
};

}
