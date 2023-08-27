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

#include <GL/glew.h>
#include "../../../include/Global.h"

namespace aprn::vis {

/** OpenGL function wrapper for error checking in debug mode. */
#ifdef DEBUG_MODE
   #define GLCall(function) { GLClearErrors(); function; GLLogCall(#function, __FILE__, __LINE__); }
#else
   #define GLCall(function) function
#endif

/** Clear all OpenGL errors. */
inline void
GLClearErrors() { while(glGetError() != GL_NO_ERROR) {} }

/** Throw current OpenGL error and exit program. */
inline void
GLLogCall(const char* function, const char* file, const int line)
{
  if(GLenum error_code = glGetError())
  {
    std::string error_message;
    switch(error_code)
    {
      case GL_INVALID_ENUM:                  error_message = "INVALID_ENUM";                  break;
      case GL_INVALID_VALUE:                 error_message = "INVALID_VALUE";                 break;
      case GL_INVALID_OPERATION:             error_message = "INVALID_OPERATION";             break;
      case GL_STACK_OVERFLOW:                error_message = "STACK_OVERFLOW";                break;
      case GL_STACK_UNDERFLOW:               error_message = "STACK_UNDERFLOW";               break;
      case GL_OUT_OF_MEMORY:                 error_message = "OUT_OF_MEMORY";                 break;
      case GL_INVALID_FRAMEBUFFER_OPERATION: error_message = "INVALID_FRAMEBUFFER_OPERATION"; break;
    }

    EXIT_FROM(file, line, "OpenGL error 0x0", std::hex, error_code, std::dec, " (", error_message,") from call:\n\t", function)
  }
}

inline void
GLDebugOutput(GLenum source, GLenum type, unsigned id, GLenum severity, [[maybe_unused]] GLsizei length,
              const char* message, [[maybe_unused]] const void* userParam)
{
   // Ignore insignificant error/warning codes
   if(id == OneOf(131169u, 131185u, 131218u, 131204u)) return;

   std::cout << "---------------" << std::endl;
   std::cout << "Debug message (" << id << "): " << message << std::endl;

   // Write error source
   switch(source)
   {
      case GL_DEBUG_SOURCE_API:               std::cout << "Source: API"; break;
      case GL_DEBUG_SOURCE_WINDOW_SYSTEM:     std::cout << "Source: Window System"; break;
      case GL_DEBUG_SOURCE_SHADER_COMPILER:   std::cout << "Source: Shader Compiler"; break;
      case GL_DEBUG_SOURCE_THIRD_PARTY:       std::cout << "Source: Third Party"; break;
      case GL_DEBUG_SOURCE_APPLICATION:       std::cout << "Source: Application"; break;
      case GL_DEBUG_SOURCE_OTHER:             std::cout << "Source: Other"; break;
   }
   std::cout << std::endl;

   // Write error type
   switch(type)
   {
      case GL_DEBUG_TYPE_ERROR:               std::cout << "Type: Error"; break;
      case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: std::cout << "Type: Deprecated Behaviour"; break;
      case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:  std::cout << "Type: Undefined Behaviour"; break;
      case GL_DEBUG_TYPE_PORTABILITY:         std::cout << "Type: Portability"; break;
      case GL_DEBUG_TYPE_PERFORMANCE:         std::cout << "Type: Performance"; break;
      case GL_DEBUG_TYPE_MARKER:              std::cout << "Type: Marker"; break;
      case GL_DEBUG_TYPE_PUSH_GROUP:          std::cout << "Type: Push Group"; break;
      case GL_DEBUG_TYPE_POP_GROUP:           std::cout << "Type: Pop Group"; break;
      case GL_DEBUG_TYPE_OTHER:               std::cout << "Type: Other"; break;
   }
   std::cout << std::endl;

   // Write error severity
   switch(severity)
   {
      case GL_DEBUG_SEVERITY_HIGH:            std::cout << "Severity: high"; break;
      case GL_DEBUG_SEVERITY_MEDIUM:          std::cout << "Severity: medium"; break;
      case GL_DEBUG_SEVERITY_LOW:             std::cout << "Severity: low"; break;
      case GL_DEBUG_SEVERITY_NOTIFICATION:    std::cout << "Severity: notification"; break;
   }
   std::cout << std::endl;
}

}
