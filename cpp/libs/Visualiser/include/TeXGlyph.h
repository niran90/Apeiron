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
#include "Colour.h"
#include "Model.h"
#include "GlyphSheet.h"
#include "TeXObject.h"
#include "TeXSpacer.h"

#include <optional>

namespace aprn::vis {

class TeXGlyph final : public TeXObject,
                       public Model
{
 public:
   TeXGlyph() = default;

   explicit TeXGlyph(char tex_char);

   explicit TeXGlyph(const std::string& tex_str);

   TeXGlyph& Set(char tex_char);

   TeXGlyph& Set(const std::string& tex_str);

   TeXGlyph& SetColour(const SVectorR4& rgba_colour) override;

   TeXGlyph& SetColour(const Colour& colour) override;

   TeXGlyph& SetItalic(bool is_italic) override;

   TeXGlyph& SetBold(bool is_bold) override;

   inline void DoNotRender() { Render_ = false; }

   inline bool Rendered() const { return Render_; }

 private:
   void InitTeXObject(GlyphSheet::IndexT& index_offset) override;

   void ComputeDimensions(const GlyphSheet& glyph_sheet, UChar font_size, const SVectorR3& texbox_anchor, const SVectorR2& texbox_dims,
                          TeXSpacer& spacer) override;

   void OffsetAlongZ(Real& cumu_offset) override;

 private:

   GlyphSheet::IndexT    Index_{MaxInt<GlyphSheet::IndexT>};
   Option<Colour>        Colour_;
   Option<bool>          Italic_;
   Option<bool>          Bold_;
   bool                  Render_{true};
   bool                  AddSpacer_{false};
   bool                  Init_{false};
   constexpr static Real OffsetZ_{1.0e-6};
};

}