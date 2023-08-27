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
#include "GlyphSheet.h"

#include <string>

namespace aprn::vis {

struct TeXSpacer
{
   inline Real Offset(const Real glyph_anchor_x, const bool spacer_after)
   {
      if(!Enabled_) return Zero;

      if(glyph_anchor_x < Threshold_) Offset_ = Zero;
      Threshold_ = glyph_anchor_x;

      const auto current_offset = -Offset_;
      if(spacer_after) Offset_ += Size_ * GlyphSheet::PointSize();
      return current_offset;
   }

   inline static std::string Text() { return Enabled_ ? "\\hspace{" + std::to_string(Size_) + "pt}" : ""; }

   inline static bool Required(const std::string_view& tex_str) { return Enabled_ && tex_str == OneOf("s", "t"); }

   inline static void Disable() { Enabled_ = false; }

 private:
   Real                 Offset_{};
   Real                 Threshold_{-LowestFloat<>};
   constexpr static int Size_{2}; // In LaTeX points (pt)
   inline static bool   Enabled_{true};
};

}
