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

#include "../include/Global.h"
#include "../../DataContainer/include/Array.h"
#include "../../LinearAlgebra/include/Vector.h"
#include "GLDebug.h"
#include "GLTypes.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

namespace aprn::vis {

class Window
{
 public:
   Window(GLint width, GLint height);

   ~Window();

   void Open();

   void Open(GLint width, GLint height);

   void InitOpenGL();

   bool isOpen() const;

   void Close();

   void Terminate();

   bool ViewportModified();

   void ResetViewport() const;

   void SetTitle(const std::string& title, const bool append = false);

   void SwapBuffers();

   void InitTime() const;

   void ComputeDeltaTime();

   void ComputeFrameRate();

   Real DeltaTime() const;

   Real CurrentTime() const;

   GLfloat ViewportAspectRatio() const;

   SVectorR2 CursorDisplacement();

   SVectorR2 WheelDisplacement();

 private:
   friend class Visualiser;

   Pair<GLint>
   ViewportDimensions() const;

   void CreateCallBacks() const;

   static void
   HandleKeys(GLFWwindow* p_window, const GLint key, const GLint code, const GLint action, const GLint mode);

   static void
   HandleMousePosition(GLFWwindow* p_window, const GLdouble x_coord, const GLdouble y_coord);

   static void
   HandleMouseWheel(GLFWwindow* p_window, const GLdouble x_offset, const GLdouble y_offset);

   constexpr static UInt nKeys {1024};

   std::string     Title_;
   Real            CurrentTime_{};
   Real            PreviousTime_{};
   Real            PreviousFpsTime_{};
   Real            DeltaTime_{};
   SVectorR2       PreviousCursorPosition_;
   SVectorR2       CursorDisplacement_;
   SVectorR2       WheelDisplacement_;
   SVector2<GLint> WindowDimensions_;
   SVector2<GLint> ViewportDimensions_;
   SArrayB<nKeys>  Keys_;
   GLFWwindow*     GLFWWindow_;
   UInt8           FrameCounter_;
   bool            FirstCursorMotion_;
};

inline void ClearFrameBuffer()
{
   GLCall(glClearColor(0.0f, 0.0f, 0.0f, 1.0f))
   GLCall(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT))
   GLCall(glEnable(GL_DEPTH_TEST))
}

}
