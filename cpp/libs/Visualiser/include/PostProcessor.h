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

#include "Buffers.h"
#include "Model.h"
#include "Shader.h"
#include "Texture.h"

namespace aprn::vis {

enum class PostProcessType
{
   HDR,
   Blur,
   MotionBlur,
   Fade,
   Vignette
};

struct FrameImage
{
   template<class type> using UMap = std::unordered_map<std::string, type>;
   FrameBuffer   FBO;
   RenderBuffer  RBO;
   UMap<Texture> Textures;
};

class PostProcessor
{
 public:
   void Init(Pair<GLint> dimensions);

   void InitMultiSampledBuffers();

   void InitResolvedBuffers();

   void InitBlurBuffers();

   void StartWrite() const;

   void StopWrite() const;

   void ResolveFrameBuffer() const;

   void Render();

 private:
   friend class Visualiser;

   template<class type> using UMap = std::unordered_map<std::string, type>;
   FrameImage                MultiSampledImage_;
   FrameImage                ResolvedImage_;
   UMap<FrameImage>          BlurBuffers_;
   UMap<Shader>              Shaders_;
   Model                     ScreenQuad_;
   UInt                      Width_;
   UInt                      Height_;
   bool                      Init_{};
   constexpr static UInt     SampleCount_{32}; // For MSAA
   inline static std::string Default_ = "Default";
};

}
