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

#include "../../../include/Global.h"
#include "../include/Window.h"
#include "LinearAlgebra/include/Vector.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Public Interface
***************************************************************************************************************************************************************/
Window::Window(GLint width, GLint height)
   : WindowDimensions_{width, height} {}

Window::~Window() { Terminate(); }

void
Window::Open() { Open(WindowDimensions_.x(), WindowDimensions_.y()); }

void
Window::Open(const GLint width, const GLint height)
{
   if(!glfwInit())
   {
      glfwTerminate();
      EXIT("Failed to Initialise GLFW.")
   }

   // Window properties
   glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
   glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
   glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // Enforces backward incompatibility
   glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
//   glfwWindowHint(GLFW_SAMPLES, 24); // Anti-aliasing multi-sample count
#ifdef DEBUG_MODE
   glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
#endif

   // Anti-aliasing properties
   glfwWindowHint(GLFW_SAMPLES, 32);
   GLCall(glEnable(GL_MULTISAMPLE))

   // Create a window and its OpenGL context.
   WindowDimensions_ = {width, height };
   GLFWWindow_ = glfwCreateWindow(width, height, "Apeiron", nullptr, nullptr);
   if(!GLFWWindow_)
   {
      glfwTerminate();
      EXIT("Could not create an OpenGL window.")
   }

   // Set context for GLEW to use
   glfwMakeContextCurrent(GLFWWindow_);
   glfwSwapInterval(1);

   // Set viewport dimensions
   std::tie(ViewportDimensions_[0], ViewportDimensions_[1]) = ViewportDimensions();

   // Handle key mouse inputs
   CreateCallBacks();
   glfwSetInputMode(GLFWWindow_, GLFW_CURSOR, GLFW_CURSOR_NORMAL); // Show cursor in window.

   // Allow modern extension features
   glewExperimental = GL_TRUE;

   // Initialise GLEW
   if(glewInit() != GLEW_OK)
   {
      Terminate();
      EXIT("Failed to Initialise GLEW.")
   }
   else Print("\nRunning OpenGL Version:", glGetString(GL_VERSION));

   glfwSetWindowUserPointer(GLFWWindow_, this);
}

void
Window::InitOpenGL()
{
   // Initialise OpenGL debug output
#ifdef DEBUG_MODE
   int flags;
   GLCall(glGetIntegerv(GL_CONTEXT_FLAGS, &flags))
   if(flags & GL_CONTEXT_FLAG_DEBUG_BIT)
   {
      GLCall(glEnable(GL_DEBUG_OUTPUT))
      GLCall(glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS))
      GLCall(glDebugMessageCallback(GLDebugOutput, nullptr))
      GLCall(glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE))
   }
#endif

   GLCall(glEnable(GL_DEPTH_TEST))
   GLCall(glEnable(GL_MULTISAMPLE))
//   GLCall(glEnable(GL_CULL_FACE))
//   GLCall(glCullFace(GL_FRONT))
//   GLCall(glFrontFace(GL_CCW))
   GLCall(glEnable(GL_BLEND))
   GLCall(glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA))
//   GLCall(glPolygonMode(GL_FRONT_AND_BACK, GL_LINE))
}

bool
Window::isOpen() const { return !glfwWindowShouldClose(GLFWWindow_); }

void
Window::Close() { glfwSetWindowShouldClose(GLFWWindow_, GL_TRUE); }

void
Window::Terminate()
{
   glfwDestroyWindow(GLFWWindow_);
   glfwTerminate();
}

bool
Window::ViewportModified()
{
   DEBUG_ASSERT(isOpen(), "The OpenGL window is not open.")

   // Adjust viewport, if necessary
   auto [width, height] = ViewportDimensions();
   if(width != ViewportDimensions_[0] || height != ViewportDimensions_[1])
   {
      ViewportDimensions_[0] = width;
      ViewportDimensions_[1] = height;
      GLCall(glViewport(0, 0, width, height));
      return true;
   }
   else return false;
}

void
Window::ResetViewport() const { GLCall(glViewport(0, 0, ViewportDimensions_[0], ViewportDimensions_[1])); }

void
Window::SetTitle(const std::string& title, const bool append)
{
   if(append) glfwSetWindowTitle(GLFWWindow_, (Title_ + title).c_str());
   else
   {
      Title_ = title;
      glfwSetWindowTitle(GLFWWindow_, Title_.c_str());
   }
}

void
Window::SwapBuffers() { glfwSwapBuffers(GLFWWindow_); }

void
Window::InitTime() const { glfwSetTime(Zero); }

void
Window::ComputeDeltaTime()
{
   CurrentTime_  = glfwGetTime();
   DeltaTime_    = CurrentTime_ - PreviousTime_;
   PreviousTime_ = CurrentTime_;
}

void
Window::ComputeFrameRate()
{
   CurrentTime_  = glfwGetTime();
   const auto delta_time = CurrentTime_ - PreviousFpsTime_;
   ++FrameCounter_;

   // Re-compute FPS.
   if(delta_time > Tenth)
   {
      const auto fps = static_cast<Real>(FrameCounter_) / delta_time;
      const auto frame_duration = 1.0e3 / fps; // in milliseconds
      const std::string title_suffix = "  |  " + ToString(fps, 2) + " fps  |  " + ToString(frame_duration, 2) + " ms";
      SetTitle(title_suffix, true);

      PreviousFpsTime_ = CurrentTime_;
      FrameCounter_ = 0;
   }
}

Real
Window::CurrentTime() const { return CurrentTime_; }

Real
Window::DeltaTime() const { return DeltaTime_; }

GLfloat
Window::ViewportAspectRatio() const { return static_cast<GLfloat>(ViewportDimensions_[0]) / static_cast<GLfloat>(ViewportDimensions_[1]); }

SVectorR2
Window::CursorDisplacement()
{
   const auto x_disp = CursorDisplacement_[0];
   const auto y_disp = CursorDisplacement_[1];
   CursorDisplacement_[0] = Zero;
   CursorDisplacement_[1] = Zero;
   return { x_disp, y_disp };
}

SVectorR2
Window::WheelDisplacement()
{
   const auto x_disp = WheelDisplacement_[0];
   const auto y_disp = WheelDisplacement_[1];
   WheelDisplacement_[0] = Zero;
   WheelDisplacement_[1] = Zero;
   return { x_disp, y_disp };
}

/***************************************************************************************************************************************************************
* Private Interface
***************************************************************************************************************************************************************/
Pair<GLint>
Window::ViewportDimensions() const
{
   GLint width, height;
   glfwGetFramebufferSize(GLFWWindow_, &width, &height);
   return { width, height };
}

void
Window::CreateCallBacks() const
{
   glfwSetKeyCallback(GLFWWindow_, HandleKeys);
   glfwSetCursorPosCallback(GLFWWindow_, HandleMousePosition);
   glfwSetScrollCallback(GLFWWindow_, HandleMouseWheel);
}

void
Window::HandleKeys(GLFWwindow* p_window, const GLint key, [[maybe_unused]] const GLint code, const GLint action, [[maybe_unused]] const GLint mode)
{
   // Get pointer to the Window object which contains p_window.
   Window* window = static_cast<Window*>(glfwGetWindowUserPointer(p_window));

   // Check for pressing of the ESC key.
   if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(p_window, GL_TRUE);

   // Check for press/release of all other keys.
   if(isBounded(key, static_cast<GLint>(0), static_cast<GLint>(nKeys)))
   {
      if     (action == GLFW_PRESS)   window->Keys_[key] = true;
      else if(action == GLFW_RELEASE) window->Keys_[key] = false;
   }
}

void
Window::HandleMousePosition(GLFWwindow* p_window, const GLdouble x_coord, const GLdouble y_coord)
{
   Window* window = static_cast<Window*>(glfwGetWindowUserPointer(p_window));

   if(window->FirstCursorMotion_)
   {
      window->PreviousCursorPosition_ = {x_coord, y_coord };
      window->FirstCursorMotion_ = false;
   }

   window->CursorDisplacement_     = {x_coord - window->PreviousCursorPosition_[0], y_coord - window->PreviousCursorPosition_[1] };
   window->PreviousCursorPosition_ = {x_coord, y_coord };
}

void
Window::HandleMouseWheel(GLFWwindow* p_window, [[maybe_unused]] const GLdouble x_offset, const GLdouble y_offset)
{
   Window* window = static_cast<Window*>(glfwGetWindowUserPointer(p_window));
   window->WheelDisplacement_ = {0.0, y_offset };
}

}

