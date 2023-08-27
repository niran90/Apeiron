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
#include "LinearAlgebra/include/Vector.h"
#include "GLDebug.h"
#include "GLTypes.h"

#include <memory>
#include <optional>
#include <GL/glew.h>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Texture Types
***************************************************************************************************************************************************************/
enum class TextureType
{
   /** Texture file types */
   Diffuse,
   Normal,
   Displacement,
   Roughness,
   AmbientOcclusion,

   /** Frame buffer object types */
   DirectionalDepth,
   PointDepth
};

/***************************************************************************************************************************************************************
* Texture Class
***************************************************************************************************************************************************************/
class Texture
{
 public:
   Texture(const TextureType type, const std::string& file_path);

   Texture(const TextureType type, const bool is_fbo_attachment);

   Texture(const Texture& texture) = delete;

   Texture(Texture&& texture) noexcept;

   ~Texture();

   void Read(const std::string& file_path, const GLint wrap_type);

   void Init(GLuint width, GLuint height, GLint internal_format, GLenum format, GLenum data_type, GLint wrap_type, size_t n_samples = 1,
             const SVector4<GLfloat>& border_colour = { 1.0f, 1.0f, 1.0f, 1.0f });

   void Bind(UInt slot = 0) const;

   void Unbind() const;

   void Delete();

   Texture& operator=(const Texture& texture) = delete;

   Texture& operator=(Texture&& texture) noexcept;

   inline auto ID() const { return ID_; }

   inline auto Width() const { return Width_; }

   inline auto Height() const { return Height_; }

   inline auto MapScale() const { return MapScale_ ? MapScale_.value() : throw "Attempting to access a map scale that has not been set."; }

   inline void SetMapScale(const Real scale) { MapScale_ = scale; }

 private:
   friend class ModelGroup;
   friend class Shadow;

   GLint OpenGLType() const;

   UInt         ID_{};
   TextureType  Type_;
   UPtr<UChar>  LocalBuffer_{};
   Option<Real> MapScale_;
   int          Width_;
   int          Height_;
   int          ChannelCount_;
   size_t       SampleCount_;
   bool         MultiSampled_;
   bool         FBOAttachment_;
};

/***************************************************************************************************************************************************************
* Texture Stand-alone Functions
***************************************************************************************************************************************************************/
std::string
TextureName(const std::string& material, const std::string& item, const size_t index, const size_t resolution);

std::string
TextureTypeString(TextureType type);

TextureType
GetTextureType(const std::string& type_string);

std::string
TextureUniformString(const TextureType type);

std::string
TextureUniformString(const std::string& type_string);

std::string
TextureDirectory(const std::string& name);

Option<std::string>
TexturePath(const std::string& file_directory, TextureType type);

}
