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
#include "DataContainer/include/Array.h"
#include "Colour.h"
#include "GlyphSheet.h"
#include "TeXSpacer.h"

#include <string>

namespace aprn::vis {

class TeXObject
{
 public:
   virtual ~TeXObject() = default;

   virtual TeXObject& SetColour(const SVectorR4& rgba_colour) = 0;

   virtual TeXObject& SetColour(const Colour& colour) = 0;

   virtual TeXObject& SetItalic(bool is_italic) = 0;

   virtual TeXObject& SetBold(bool is_bold) = 0;

   virtual void ComputeDimensions(const GlyphSheet& glyph_sheet, UChar font_size, const SVectorR3& texbox_anchor, const SVectorR2& texbox_dimensions,
                                  TeXSpacer& spacer) = 0;

   virtual void InitTeXObject(GlyphSheet::IndexT& index_offset) = 0;

   virtual void OffsetAlongZ(Real& cumu_offset) = 0;

   inline void SetText(const std::string& text) { Text_ = text; }

   inline void AddText(const std::string& text) { Text_ += text; }

   inline const auto& Text() const { return Text_; }

 protected:
   std::string Text_;
};

}
