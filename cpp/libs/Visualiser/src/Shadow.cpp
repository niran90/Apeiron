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

#include "../include/Shadow.h"

namespace aprn::vis {

Shadow::Shadow(const bool is_point_source)
   : FBO_(), DepthMap_(is_point_source ? TextureType::PointDepth : TextureType::DirectionalDepth, true), PointSource_(is_point_source) {}

void Shadow::Init(const GLsizei width, const GLsizei height)
{
   DepthMap_.Init(width, height, GL_DEPTH_COMPONENT, GL_DEPTH_COMPONENT, GL_FLOAT, PointSource_ ? GL_CLAMP_TO_EDGE : GL_CLAMP_TO_BORDER);

   FBO_.Init();
   FBO_.Bind();
   PointSource_ ? FBO_.AttachTexture(GL_DEPTH_ATTACHMENT, DepthMap_.ID()) :
                  FBO_.AttachTexture2D(GL_DEPTH_ATTACHMENT, DepthMap_.ID());
   FBO_.Draw(GL_NONE);
   FBO_.Read(GL_NONE);
   FBO_.Check();
   FBO_.Unbind();
}

void Shadow::StartWrite() const
{
   FBO_.Bind();
   GLCall(glClear(GL_DEPTH_BUFFER_BIT))
}

}

