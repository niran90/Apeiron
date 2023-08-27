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
#include "../../LinearAlgebra/include/Vector.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

template<typename T>
consteval GLenum
GLTypeEnum([[maybe_unused]] const T& value = T())
{
  if constexpr     (isTypeSame<T, GLbyte>())   return GL_BYTE;
  else if constexpr(isTypeSame<T, GLubyte>())  return GL_UNSIGNED_BYTE;
  else if constexpr(isTypeSame<T, GLshort>())  return GL_SHORT;
  else if constexpr(isTypeSame<T, GLushort>()) return GL_UNSIGNED_SHORT;
  else if constexpr(isTypeSame<T, GLint>())    return GL_INT;
  else if constexpr(isTypeSame<T, GLuint>())   return GL_UNSIGNED_INT;
  else if constexpr(isTypeSame<T, GLfloat>())  return GL_FLOAT;
  else if constexpr(isTypeSame<T, GLdouble>()) return GL_DOUBLE;
  else EXIT("Unrecognised OpenGL data type: ")
}

constexpr GLuint
GLTypeSize(const GLenum gl_type)
{
  switch(gl_type)
  {
    case GL_BYTE:           return sizeof(GLbyte);
    case GL_UNSIGNED_BYTE:  return sizeof(GLubyte);
    case GL_SHORT:          return sizeof(GLshort);
    case GL_UNSIGNED_SHORT: return sizeof(GLushort);
    case GL_INT:            return sizeof(GLint);
    case GL_UNSIGNED_INT:   return sizeof(GLuint);
    case GL_FLOAT:          return sizeof(GLfloat);
    case GL_DOUBLE:         return sizeof(GLdouble);
  }
  EXIT("Unrecognised OpenGL data type.")
}

/** Convert a GLM vector to a static vector. */
template<typename T, size_t N, glm::qualifier Q = glm::defaultp>
consteval SVector<T, N>
GlmVecToSVector(const glm::vec<N, T, Q>& in)
{
  SVector<T, N> out;
  FOR(i, N) out[i] = in[i];
  return out;
}

/** Convert a static vector to a GLM vector. */
template<typename T, size_t N, glm::qualifier Q = glm::defaultp>
constexpr glm::vec<N, T, Q>
SVectorToGlmVec(const SVector<T, N>& in)
{
  glm::vec<N, T, Q> out;
  FOR(i, N) out[i] = in[i];
  return out;
}

}
