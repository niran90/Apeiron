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

   void CursorControl(const SVectorR2& cursor_displacement);

   void WheelControl(const SVectorR2& wheel_displacement);

   void UpdateViewMatrix();

   void UpdateProjMatrix();

   inline const glm::vec3& Position() const { return Position_; }

   inline const glm::mat4& ViewMatrix() const { return ViewMatrix_; }

   inline const glm::mat4& ProjMatrix() const { return ProjMatrix_; }

 private:
   friend class Visualiser;

   void ClipPitch() { Clip(Pitch_, -89.0f, 89.0f); }

   glm::mat4 ViewMatrix_{1.0};
   glm::mat4 ProjMatrix_{1.0};
   glm::vec3 Position_;
   glm::vec3 Front_;
   glm::vec3 Up_;
   glm::vec3 Right_;
   GLfloat   Pitch_;
   GLfloat   Yaw_;
   GLfloat   AspectRatio_;
   GLfloat   FieldOfView_;
   GLfloat   NearPlane_;
   GLfloat   FarPlane_;

   constexpr static glm::vec3 WorldUp_ {0.0f, 1.0f, 0.0f};
   constexpr static GLfloat   MotionSensitivity_   {6.0f};
   constexpr static GLfloat   RotationSensitivity_ {7.5f};
   constexpr static GLfloat   CursorSensitivity_   {0.4f};
   constexpr static GLfloat   WheelSensitivity_    {1.0f};
};

}
