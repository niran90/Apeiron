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

#include "../include/Texture.h"
#include "FileManager/include/FileSystem.h"

#include <GLFW/glfw3.h>
#include "../../../external/stb_image/stb_image.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Texture Class Implementation
***************************************************************************************************************************************************************/
Texture::Texture(const TextureType type, const std::string& file_path)
   : Texture(type, false)
{
   Read(file_path, GL_CLAMP_TO_EDGE);
}

Texture::Texture(const TextureType type, const bool is_fbo_attachment)
  : ID_(0), Type_(type), LocalBuffer_(nullptr), Width_(0), Height_(0), ChannelCount_(0), FBOAttachment_(is_fbo_attachment) {}

Texture::Texture(Texture&& texture) noexcept { *this = std::move(texture); }

Texture::~Texture() { Delete(); }

void
Texture::Read(const std::string& file_path, const GLint wrap_type)
{
   // Load texture image with stbi.
   int width, height;
   stbi_set_flip_vertically_on_load(true);
   LocalBuffer_.reset(stbi_load(file_path.c_str(), &width, &height, &ChannelCount_, 0));
   ASSERT(LocalBuffer_, "Could not load file \"", file_path, "\" to texture.")

   // Initialise texture.
   std::pair<GLint, GLenum> format = ChannelCount_ == 2 ? std::make_pair(GL_SRGB, GL_RG) :
                                     ChannelCount_ == 3 ? std::make_pair(GL_SRGB, GL_RGB) :
                                     ChannelCount_ == 4 ? std::make_pair(GL_SRGB_ALPHA, GL_RGBA) :
                                     throw "Currently only 2, 3, or 4 bits per pixel are supported. Got " + ToString(ChannelCount_) + ".";
   Init(width, height, format.first, format.second, GL_UNSIGNED_BYTE, wrap_type);

   // Free stbi-loaded image buffer.
   stbi_image_free(LocalBuffer_.get());
}

void
Texture::Init(const GLuint width, const GLuint height, const GLint internal_format, const GLenum format, const GLenum data_type,
              const GLint wrap_type, const size_t n_samples, const SVector4<GLfloat>& border_colour)
{
   ASSERT(glfwGetCurrentContext(), "Cannot intialise a texture without an OpenGL context.")
   ASSERT(n_samples > 0, "The sample count for each texture must be at least one.")

   GLCall(glGenTextures(1, &ID_))

   Width_        = width;
   Height_       = height;
   SampleCount_  = n_samples;
   MultiSampled_ = n_samples > 1;

   Bind();

   // Common settings
   const auto gl_type = OpenGLType();
   if(!MultiSampled_)
   {
      GLCall(glTexParameteri(gl_type, GL_TEXTURE_MIN_FILTER, GL_LINEAR))
      GLCall(glTexParameteri(gl_type, GL_TEXTURE_MAG_FILTER, GL_LINEAR))
   }

   // OpenGL type-specific settings
   if(gl_type == OneOf(GL_TEXTURE_2D, GL_TEXTURE_2D_MULTISAMPLE))
   {
     if(wrap_type == GL_CLAMP_TO_BORDER) GLCall(glTexParameterfv(gl_type, GL_TEXTURE_BORDER_COLOR, border_colour.data()))

     if(gl_type == GL_TEXTURE_2D)
     {
        GLCall(glTexParameteri(gl_type, GL_TEXTURE_WRAP_S, wrap_type))
        GLCall(glTexParameteri(gl_type, GL_TEXTURE_WRAP_T, wrap_type))
        GLCall(glTexImage2D(gl_type, 0, internal_format, Width_, Height_, 0, format, data_type, LocalBuffer_.get()))
     }
     else GLCall(glTexImage2DMultisample(gl_type, SampleCount_, internal_format, Width_, Height_, GL_TRUE))
   }
   else if(gl_type == GL_TEXTURE_CUBE_MAP)
   {
     ASSERT(wrap_type == GL_CLAMP_TO_EDGE, "Only GL_CLAMP_TO_EDGE is currently supported for cube maps.")

     GLCall(glTexParameteri(gl_type, GL_TEXTURE_WRAP_S, wrap_type))
     GLCall(glTexParameteri(gl_type, GL_TEXTURE_WRAP_T, wrap_type))
     GLCall(glTexParameteri(gl_type, GL_TEXTURE_WRAP_R, wrap_type))

     FOR(i, 6) GLCall(glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, internal_format, Width_, Height_, 0, format, data_type, nullptr))
   }

   if(!FBOAttachment_) GLCall(glGenerateMipmap(gl_type))

   Unbind();
}

void
Texture::Bind(const UInt slot) const
{
   GLCall(glActiveTexture(GL_TEXTURE0 + slot))
   GLCall(glBindTexture(OpenGLType(), ID_))
}

void
Texture::Unbind() const { GLCall(glBindTexture(OpenGLType(), 0)) }

void
Texture::Delete()
{
   GLCall(glDeleteTextures(1, &ID_))
   ID_ = 0;
}

Texture&
Texture::operator=(Texture&& texture) noexcept
{
   // Move all memebers.
   ID_            = std::move(texture.ID_);
   Type_          = std::move(texture.Type_);
   LocalBuffer_   = std::move(texture.LocalBuffer_);
   Width_         = std::move(texture.Width_);
   Height_        = std::move(texture.Height_);
   ChannelCount_  = std::move(texture.ChannelCount_);
   MapScale_      = std::move(texture.MapScale_);
   SampleCount_   = std::move(texture.SampleCount_);
   MultiSampled_  = std::move(texture.MultiSampled_);
   FBOAttachment_ = std::move(texture.FBOAttachment_);

   // Nullify ID of moved texture.
   texture.ID_ = 0;

   return *this;
}

GLint
Texture::OpenGLType() const
{
   typedef TextureType TT;
   if(Type_ == OneOf(TT::Diffuse, TT::Normal, TT::Displacement, TT::DirectionalDepth)) return MultiSampled_ ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D;
   else if(Type_ == TT::PointDepth) return GL_TEXTURE_CUBE_MAP;
   else EXIT("Unrecongnised texture type.")
}

/***************************************************************************************************************************************************************
* Texture Stand-alone Functions
***************************************************************************************************************************************************************/
std::string
TextureName(const std::string& material, const std::string& item, const size_t index, const size_t resolution)
{
   ASSERT(resolution == OneOf(1U, 2U, 4U), "Only 1K, 2K, and 4K resolutions are recognised.")
   return material + "/" + (item.empty() ? "" : item + "/") + ToString(index) + "/" + ToString(resolution) + "K";
}

std::string
TextureTypeString(const TextureType type)
{
   switch(type)
   {
      case TextureType::Diffuse:          return "Diffuse";
      case TextureType::Normal:           return "Normal";
      case TextureType::Displacement:     return "Displacement";
      case TextureType::Roughness:        return "Roughness";
      case TextureType::AmbientOcclusion: return "AmbientOcclusion";
      case TextureType::DirectionalDepth: return "DirectionalDepth";
      case TextureType::PointDepth:       return "PointDepth";
      default: EXIT("Unrecognised texture type.")
   }
}

TextureType
GetTextureType(const std::string& type_string)
{
   if     (type_string == "Diffuse")          return TextureType::Diffuse;
   else if(type_string == "Normal")           return TextureType::Normal;
   else if(type_string == "Displacement")     return TextureType::Displacement;
   else if(type_string == "Roughness")        return TextureType::Roughness;
   else if(type_string == "AmbientOcclusion") return TextureType::AmbientOcclusion;
   else if(type_string == "DirectionalDepth") return TextureType::DirectionalDepth;
   else if(type_string == "PointDepth")       return TextureType::PointDepth;
   else EXIT("Unrecognised texture type.")
}

std::string
TextureUniformString(const TextureType type)
{
   switch(type)
   {
      case TextureType::Diffuse:      return "diffuse_map";
      case TextureType::Normal:       return "normal_map";
      case TextureType::Displacement: return "displacement_map";
      default: EXIT("The given texture type does not have an associate uniform string.")
   }
}

std::string
TextureUniformString(const std::string& type_string) { return TextureUniformString(GetTextureType(type_string)); }

std::string
TextureDirectory(const std::string& name) { return "libs/Visualiser/resources/textures/" + name + "/"; }

Option<std::string>
TexturePath(const std::string& file_directory, const TextureType type)
{
   const auto& texture_str = TextureTypeString(type);
   std::string file_path = file_directory + texture_str;

   if(EnumToInt(type) <= EnumToInt(TextureType::AmbientOcclusion))
   {
      // Try finding a .png file first, failing which, find a .jpg file
      for(const auto& ext : {".png", ".jpg"})
         if(flmgr::FileExists(file_path + ext))
         {
            file_path.append(ext);
            return std::make_optional(file_path);
         }
      return std::nullopt;
   }
   else EXIT("The texture type ", texture_str, " cannot be read from file.")
}

}

