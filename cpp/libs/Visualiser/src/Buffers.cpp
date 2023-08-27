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

#include "../include/Buffers.h"
#include <GLFW/glfw3.h>
#include <cstring>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Buffer Abstract Base Class
***************************************************************************************************************************************************************/
namespace detail {

template<BufferType T>
Buffer<T>::Buffer(Buffer<T>&& buffer) noexcept { *this = std::move(buffer); }

template<BufferType T>
Buffer<T>::~Buffer() { Delete(); }

template<BufferType T>
void
Buffer<T>::Init()
{
   ASSERT(glfwGetCurrentContext(), "Cannot intialise buffer without an OpenGL context.")

   if constexpr(T == OneOf(BT::VBO, BT::EBO, BT::SSBO)) GLCall(glGenBuffers(1, &ID_))
   else if(T == BT::VAO)                                GLCall(glGenVertexArrays(1, &ID_))
   else if(T == BT::FBO)                                GLCall(glGenFramebuffers(1, &ID_))
   else if(T == BT::RBO)                                GLCall(glGenRenderbuffers(1, &ID_))
   else throw "Cannot initialise buffer - unrecognised buffer type.";
}

template<BufferType T>
void
Buffer<T>::Bind() const { Bind(ID_); }

template<BufferType T>
void
Buffer<T>::Unbind() const { Bind(0); }

template<BufferType T>
void
Buffer<T>::Delete()
{
   if(ID_ != 0)
   {
      if constexpr(T == OneOf(BT::VBO, BT::EBO, BT::SSBO)) GLCall(glDeleteBuffers(1, &ID_))
      else if constexpr(T == BT::VAO)                      GLCall(glDeleteVertexArrays(1, &ID_))
      else if constexpr(T == BT::FBO)                      GLCall(glDeleteFramebuffers(1, &ID_))
      else if constexpr(T == BT::RBO)                      GLCall(glDeleteRenderbuffers(1, &ID_))
      else throw "Cannot delete buffer - unrecognised buffer type.";
      ID_ = 0;
   }
}

template<BufferType T>
constexpr GLenum
Buffer<T>::Target() const
{
   if constexpr     (T == BT::VBO)  return GL_ARRAY_BUFFER;
   else if constexpr(T == BT::EBO)  return GL_ELEMENT_ARRAY_BUFFER;
   else if constexpr(T == BT::SSBO) return GL_SHADER_STORAGE_BUFFER;
   else if constexpr(T == BT::FBO)  return GL_FRAMEBUFFER;
   else if constexpr(T == BT::RBO)  return GL_RENDERBUFFER;
   else throw "No target for the given buffer type.";
}

template<BufferType T>
Buffer<T>&
Buffer<T>::operator=(Buffer<T>&& buffer) noexcept
{
   Delete();
   ID_ = buffer.ID_;
   buffer.ID_ = 0; // Note: must not call buffer.Delete() here!
   return *this;
}

template<BufferType T>
void
Buffer<T>::Bind(const GLuint id) const
{
   if constexpr(T == BT::VBO)       GLCall(glBindBuffer(Target(), id))
   else if constexpr(T == BT::EBO)  GLCall(glBindBuffer(Target(), id))
   else if constexpr(T == BT::SSBO) GLCall(glBindBuffer(Target(), id))
   else if constexpr(T == BT::FBO)  GLCall(glBindFramebuffer(Target() , id))
   else if constexpr(T == BT::RBO)  GLCall(glBindRenderbuffer(Target(), id))
   else if constexpr(T == BT::VAO)  GLCall(glBindVertexArray(id))
   else throw "Cannot bind/unbind buffer - unrecognised buffer type.";
}

template class Buffer<BufferType::VBO>;
template class Buffer<BufferType::VAO>;
template class Buffer<BufferType::EBO>;
template class Buffer<BufferType::FBO>;
template class Buffer<BufferType::RBO>;
template class Buffer<BufferType::SSBO>;

}

/***************************************************************************************************************************************************************
* Vertex Buffer Class
***************************************************************************************************************************************************************/
void
VertexBuffer::Init(const DArray<Vertex>& vertices)
{
   Buffer::Init();
   Bind();
   Load(vertices);
   Unbind();
}

void
VertexBuffer::Load(const DArray<Vertex>& vertices)
{
   VertexCount_ = vertices.size();
   GLCall(glBufferData(GL_ARRAY_BUFFER, VertexCount_ * sizeof(Vertex), vertices.data(), GL_DYNAMIC_DRAW))
}

void
VertexBuffer::Update(const DArray<Vertex>& vertices)
{
   Bind();

   // Get vertex buffer data.
   const auto target = Target();
   void *ptr = glMapBuffer(target, GL_WRITE_ONLY);

   // Now copy data into memory.
   std::memcpy(ptr, vertices.data(), VertexCount_ * sizeof(Vertex));

   // Tell OpenGL that we're done with the pointer.
   GLCall(glUnmapBuffer(target))

   Unbind();
}

/***************************************************************************************************************************************************************
* Vertex Array Class
***************************************************************************************************************************************************************/
void
VertexArray::AddBuffer(const VertexBuffer& vertex_buffer, const VertexAttributeLayout& vertex_layout)
{
   vertex_buffer.Bind();

   GLuint offset(0);
   FOR(i, vertex_layout.Attributes.size())
   {
      const auto& element = vertex_layout.Attributes[i];
      GLCall(glVertexAttribPointer(i, element.nComponents, element.GLType, element.isNormalised, vertex_layout.Stride, reinterpret_cast<void*>(offset)))
      GLCall(glEnableVertexAttribArray(i))

      offset += element.nComponents * GLTypeSize(element.GLType);
   }

   vertex_buffer.Unbind();
}

/***************************************************************************************************************************************************************
* Index Buffer Class
***************************************************************************************************************************************************************/
void
IndexBuffer::Init(const DArray<GLuint>& indices)
{
   Buffer::Init();
   Bind();
   Load(indices);
   Unbind();
}

void
IndexBuffer::Load(const DArray<GLuint>& indices)
{
   IndexCount_ = indices.size();
   GLCall(glBufferData(GL_ELEMENT_ARRAY_BUFFER, IndexCount_ * sizeof(GLuint), indices.data(), GL_STATIC_DRAW))
}

/***************************************************************************************************************************************************************
* Frame Buffer Class
***************************************************************************************************************************************************************/
void
FrameBuffer::Init(bool multi_sampled)
{
   MultiSampled_ = multi_sampled;
   Buffer::Init();
}

void
FrameBuffer::AttachTexture(const GLenum attachement, const GLuint texture_id) const
{
   GLCall(glFramebufferTexture(GL_FRAMEBUFFER, attachement, texture_id, 0))
   Check();
}

void
FrameBuffer::AttachTexture2D(const GLenum attachement, const GLuint texture_id) const
{
   GLCall(glFramebufferTexture2D(GL_FRAMEBUFFER, attachement, MultiSampled_ ? GL_TEXTURE_2D_MULTISAMPLE : GL_TEXTURE_2D, texture_id, 0))
   Check();
}

void
FrameBuffer::AttachRenderBuffer(const GLenum attachement, const GLuint rbo_id) const
{
   GLCall(glFramebufferRenderbuffer(GL_FRAMEBUFFER, attachement, GL_RENDERBUFFER, rbo_id))
   Check();
}

void
FrameBuffer::Check() const
{
   auto fb_status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
   if(fb_status != GL_FRAMEBUFFER_COMPLETE)
   {
      std::string prefix = "Could not intialise frame buffer. Reason: ";
      switch(fb_status) 
      {
         case GL_FRAMEBUFFER_UNDEFINED:                     EXIT(prefix, "The specified frame buffer is the default one, which does not yet exist.")
         case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         EXIT(prefix, "Some of the attachment points are incomplete.")
         case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: EXIT(prefix, "No image attachment found.")
         case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:        EXIT(prefix, "Incomplete draw buffer.")
         case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:        EXIT(prefix, "Incomplete read buffer.")
         case GL_FRAMEBUFFER_UNSUPPORTED:                   EXIT(prefix, "Some attached images have incompatible internal formats.")
         case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:        EXIT(prefix, "Incompatible multi-sampling for some attachments.")
         case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS:      EXIT(prefix, "Incomplete layer targets.")
      }
   }
}

void
FrameBuffer::Draw(const GLenum attachment) const { GLCall(glNamedFramebufferDrawBuffer(ID_, attachment)) }

void
FrameBuffer::Draw(const DArray<GLenum>& attachments) const { GLCall(glNamedFramebufferDrawBuffers(ID_, attachments.size(), attachments.data())) }

void
FrameBuffer::Read(const GLenum attachment) const { GLCall(glNamedFramebufferReadBuffer(ID_, attachment)) }

/***************************************************************************************************************************************************************
* Animate Buffer Class
***************************************************************************************************************************************************************/
void
RenderBuffer::Init(const size_t n_samples)
{
   ASSERT(n_samples > 0, "The sample count for each render buffer must be at least one.")
   SampleCount_  = n_samples;
   MultiSampled_ = n_samples > 1;
   Buffer::Init();
}

void
RenderBuffer::Allocate(const GLenum format, const GLsizei width, const GLsizei height)
{
   Bind();
   if(MultiSampled_) GLCall(glRenderbufferStorageMultisample(GL_RENDERBUFFER, SampleCount_, format, width, height))
   else              GLCall(glRenderbufferStorage(GL_RENDERBUFFER, format, width, height))
   Unbind();
}

}