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

#include "../include/PostProcessor.h"
#include "../include/ObjectFactory.h"
#include "../include/Window.h"

namespace aprn::vis {

void
PostProcessor::Init(Pair<GLint> dimensions)
{
   // Set member variables.
   std::tie(Width_, Height_) = dimensions;

   // Load required post-processing shaders.
   for(std::string shader : { "HDR", "Blur", "Blend" }) Shaders_.emplace(shader, Shader::Directory + shader + ".glsl");

   // Initialise required post-processing buffers.
   InitMultiSampledBuffers();
   InitResolvedBuffers();
   InitBlurBuffers();

   // Create rectangular screen-filling quad (a square in normalised device coordinates) to draw each texture to.
   ScreenQuad_ = ObjectFactory::ScreenQuad();
   ScreenQuad_->Init();

   Init_ = true;
}

void
PostProcessor::InitMultiSampledBuffers()
{
   auto& fbo      = MultiSampledImage_.FBO;
   auto& rbo      = MultiSampledImage_.RBO;
   auto& textures = MultiSampledImage_.Textures;

   // Initialise and bind an MSAA frame buffer.
   fbo.Init(true);
   fbo.Bind();

   // Add MSAA HDR buffer and bright buffer textures.
   textures.emplace("HDR"   , Texture(TextureType::Diffuse, true));
   textures.emplace("Bright", Texture(TextureType::Diffuse, true));

   // Attach MSAA textures as the frame buffer's colour buffers.
   size_t i = 0;
   DArray<GLenum> attachments;
   for (std::string name : { "HDR", "Bright" }) // Note: must be attached to the frame buffer in this order.
   {
      auto& texture = textures.at(name);
      attachments.push_back(GL_COLOR_ATTACHMENT0 + i++);

      texture.Init(Width_, Height_, GL_RGBA16F, GL_RGBA, GL_FLOAT, GL_CLAMP_TO_EDGE, SampleCount_);
      fbo.AttachTexture2D(attachments.back(), texture.ID());
   }
   fbo.Draw(attachments);

   // Attach a render buffer as the frame buffer's depth/stencil buffer.
   rbo.Init(SampleCount_);
   rbo.Allocate(GL_DEPTH24_STENCIL8, Width_, Height_);
   fbo.AttachRenderBuffer(GL_DEPTH_STENCIL_ATTACHMENT, rbo.ID());

   fbo.Unbind();
}

void
PostProcessor::InitResolvedBuffers()
{
   auto& fbo      = ResolvedImage_.FBO;
   auto& textures = ResolvedImage_.Textures;

   // Initialise and bind a non-MSAA frame buffer.
   fbo.Init(false);
   fbo.Bind();

   // Add non-MSAA HDR buffer and bright buffer textures.
   textures.emplace("HDR"   , Texture(TextureType::Diffuse, true));
   textures.emplace("Bright", Texture(TextureType::Diffuse, true));

   // Attach textures as frame buffer's colour buffers.
   size_t i = 0;
   DArray<GLenum> attachments;
   for (std::string name : { "HDR", "Bright" }) // Note: must be attached to the frame buffer in this order.
   {
      auto& texture = textures.at(name);
      attachments.push_back(GL_COLOR_ATTACHMENT0 + i++);

      texture.Init(Width_, Height_, GL_RGBA16F, GL_RGBA, GL_FLOAT, GL_CLAMP_TO_EDGE);
      fbo.AttachTexture2D(attachments.back(), texture.ID());
   }
   fbo.Draw(attachments);

   fbo.Unbind();
}

void
PostProcessor::InitBlurBuffers()
{
   // Configure ping-pong buffers
   for(auto name : { "Ping", "Pong" })
   {
      auto& buffer   = BlurBuffers_[name]; // Note: emplacement is intended here.
      auto& fbo      = buffer.FBO;
      auto& textures = buffer.Textures;

      // Initialise and bind a non-MSAA frame buffer.
      fbo.Init(false);
      fbo.Bind();

      // Attach a single colour buffer texture.
      textures.emplace(Default_, Texture(TextureType::Diffuse, true));
      auto& texture = textures.at(Default_);
      texture.Init(Width_, Height_, GL_RGBA16F, GL_RGBA, GL_FLOAT, GL_CLAMP_TO_EDGE);
      fbo.AttachTexture2D(GL_COLOR_ATTACHMENT0, texture.ID());
      fbo.Draw(GL_COLOR_ATTACHMENT0);

      fbo.Unbind();
   }
}

void
PostProcessor::StartWrite() const
{
   DEBUG_ASSERT(Init_, "The postprocessor has not yet been initialised.")
   MultiSampledImage_.FBO.Bind();
   ClearFrameBuffer();
}

void
PostProcessor::StopWrite() const
{
   DEBUG_ASSERT(Init_, "The postprocessor has not yet been initialised.")
   ResolveFrameBuffer();
   MultiSampledImage_.FBO.Unbind();
}

void
PostProcessor::ResolveFrameBuffer() const
{
   glBindFramebuffer(GL_READ_FRAMEBUFFER, MultiSampledImage_.FBO.ID());
   glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ResolvedImage_.FBO.ID());
   FOR(i, 2)
   {
      const GLenum attachment = GL_COLOR_ATTACHMENT0 + i;
      glReadBuffer(attachment);
      glDrawBuffer(attachment);
      glBlitFramebuffer(0, 0, Width_, Height_, 0, 0, Width_, Height_, GL_COLOR_BUFFER_BIT, GL_NEAREST);
   }
}

void
PostProcessor::Render()
{
   DEBUG_ASSERT(Init_, "The postprocessor has not yet been initialised.")

   GLCall(glClearColor(1.0f, 1.0f, 1.0f, 1.0f))
   GLCall(glClear(GL_COLOR_BUFFER_BIT))
   GLCall(glDisable(GL_DEPTH_TEST))

   // Apply Gaussian blur.
   DArray<FrameImage*> buffers = { &BlurBuffers_.at("Ping"),
                                   &BlurBuffers_.at("Pong") };
   const size_t n_passes = 6;
   bool horizontal = true;
   auto& blur_shader = Shaders_.at("Blur");
   blur_shader.Bind();
   for(size_t i = 0; i < n_passes; i++)
   {
      auto ping = buffers[ horizontal];
      auto pong = buffers[!horizontal];

      // Bind the FBO of the Ping buffer and the texture of the Pong buffer.
      auto& texture = i == 0 ? ResolvedImage_.Textures.at("Bright") : pong->Textures.at(Default_);
      ping->FBO.Bind();
      texture.Bind();

      blur_shader.UseTexture(texture, "u_texture", 0);
      blur_shader.SetUniform1i("u_horizontal", horizontal);
      ScreenQuad_->Render(blur_shader);

      texture.Unbind();
      ping->FBO.Unbind();

      horizontal = !horizontal;
   }
   blur_shader.Unbind();

   // Now blend the HDR and blurred bright textures.
   auto& hdr_texture  = ResolvedImage_.Textures.at("HDR");
   auto& blur_texture = buffers[!horizontal]->Textures.at(Default_);
   auto& blend_shader = Shaders_.at("Blend");
   blend_shader.Bind();

   // Load textures.
   blend_shader.UseTexture(hdr_texture , "u_hdr_texture" , 0);
   blend_shader.UseTexture(blur_texture, "u_blur_texture", 1);
   blend_shader.SetUniform1f("u_exposure", 1.0);

   ScreenQuad_->Render(blend_shader);

   hdr_texture.Unbind();
   blur_texture.Unbind();
   blend_shader.Unbind();
}

}
