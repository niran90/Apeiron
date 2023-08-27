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
#include "GLDebug.h"
#include "GLTypes.h"
#include "Buffers.h"
#include "Texture.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

class Shadow
{
 public:
   Shadow(const bool is_point_source);

   Shadow(const Shadow& shadow) = delete;

   Shadow(Shadow&& shadow) noexcept = default;

   void Init(const GLsizei width, const GLsizei height);

   void StartWrite() const;

   inline void StopWrite() const { FBO_.Unbind(); }

   inline void StartRead(const UInt slot) const { DepthMap_.Bind(slot); }

   inline const Texture& DepthMap() const { return DepthMap_; }

   Shadow& operator=(const Shadow& shadow) = delete;

   Shadow& operator=(Shadow&& shadow) noexcept = default;

 protected:
   friend class Light;

   FrameBuffer FBO_;
   Texture     DepthMap_;
   bool        PointSource_;
};

}

