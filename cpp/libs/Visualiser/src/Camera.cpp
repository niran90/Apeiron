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

#include "../include/Camera.h"

namespace aprn::vis {

Camera::Camera()
  : Camera(glm::vec3(0.0f, 0.0f, 0.0), 0.0f, 0.0f) {}

Camera::Camera(const glm::vec3& position, GLfloat pitch, GLfloat yaw)
  : Position_(position), Front_(glm::vec3(0.0f, 0.0f, -1.0f)), Pitch_(pitch), Yaw_(yaw), AspectRatio_(0.0), FieldOfView_(45.0), NearPlane_(0.0),
    FarPlane_(0.0) { SetOrientation(Position_, Pitch_, Yaw_); }

void Camera::SetOrientation(const glm::vec3& position, GLfloat pitch, GLfloat yaw)
{
   Position_ = position;

   pitch = glm::radians(pitch);
   yaw   = glm::radians(yaw);

   Front_.x = glm::cos(yaw) * glm::cos(pitch);
   Front_.y = glm::sin(pitch);
   Front_.z = glm::sin(yaw) * glm::cos(pitch);
   Front_   = glm::normalize(Front_);

   Right_ = glm::normalize(glm::cross(Front_, WorldUp_));
   Up_    = glm::normalize(glm::cross(Right_, Front_));
}

void Camera::SetViewFrustum(const GLfloat& aspect_ratio, const GLfloat& field_of_view, const GLfloat& near_plane, const GLfloat& far_plane)
{
   AspectRatio_ = aspect_ratio;
   if(field_of_view > 0.0)
   {
      FieldOfView_ = field_of_view;
      NearPlane_   = near_plane;
      FarPlane_    = far_plane;
   }
}

void Camera::KeyControl(const StaticArray<Bool, KeyCount>& keys, const GLfloat& delta_time)
{
   const auto displacement = MotionSensitivity_ * delta_time;

   // Position control
   if(keys[GLFW_KEY_W]) Position_ += displacement * Front_;
   if(keys[GLFW_KEY_S]) Position_ -= displacement * Front_;
   if(keys[GLFW_KEY_A]) Position_ -= displacement * Right_;
   if(keys[GLFW_KEY_D]) Position_ += displacement * Right_;
   if(keys[GLFW_KEY_Q]) Position_ += displacement * WorldUp_;
   if(keys[GLFW_KEY_E]) Position_ -= displacement * WorldUp_;

   // Pitch control
   if(keys[GLFW_KEY_UP])
   {
      Pitch_ += RotationSensitivity_ * displacement;
      ClipPitch();
   }
   if(keys[GLFW_KEY_DOWN])
   {
      Pitch_ -= RotationSensitivity_ * displacement;
      ClipPitch();
   }

   // Yaw control
   if(keys[GLFW_KEY_LEFT])  Yaw_ -= RotationSensitivity_ * displacement;
   if(keys[GLFW_KEY_RIGHT]) Yaw_ += RotationSensitivity_ * displacement;

   // Update camera orientation
   SetOrientation(Position_, Pitch_, Yaw_);
}

void Camera::CursorControl(const SVectorR2& cursor_displacement)
{
   // Update yaw and pitch. Ensure that the pitch is in the range [-90, 90]
   Yaw_   += CursorSensitivity_ * cursor_displacement[0];
   Pitch_ -= CursorSensitivity_ * cursor_displacement[1];
   ClipPitch();

   // Update camera orientation
   SetOrientation(Position_, Pitch_, Yaw_);
}

void Camera::WheelControl(const SVectorR2& wheel_displacement)
{
   GLfloat displacement = WheelSensitivity_ * wheel_displacement[1];
   Position_ += displacement * Front_;
}

void Camera::UpdateViewMatrix()
{
//  ViewMatrix = glm::lookAt(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f), WorldUp);
   ViewMatrix_ = glm::lookAt(Position_, Position_ + Front_, WorldUp_);
}

void Camera::UpdateProjMatrix()
{
//  ProjMatrix_ = glm::ortho(-AspectRatio, AspectRatio, -1.0f, 1.0f, -10.0f, 10.0f);
   ProjMatrix_ = glm::perspective(glm::radians(FieldOfView_), AspectRatio_, NearPlane_, FarPlane_);
}

}

