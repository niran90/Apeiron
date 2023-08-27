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
#include "ModelGroup.h"
#include "TeXObject.h"
#include "TeXGlyph.h"
#include "TeXSpacer.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeXBox Class
***************************************************************************************************************************************************************/
class TeXBox final : public TeXObject,
                     public ModelGroup
{
 public:
   TeXBox() = default;

   explicit TeXBox(const std::string& str);

   explicit TeXBox(const DArray<TeXGlyph>& tex_glyphs);

   explicit TeXBox(const DArray<TeXBox>& tex_boxes);

   explicit TeXBox(DArray<TeXGlyph>&& tex_glyphs);

   explicit TeXBox(DArray<TeXBox>&& tex_boxes);

   TeXBox& Add(const std::string& str);

   TeXBox& Add(const TeXGlyph& tex_glyph);

   TeXBox& Add(const TeXBox& tex_box);

   TeXBox& Add(const SPtr<TeXObject>& tex_object);

   TeXBox& Add(const DArray<TeXGlyph>& tex_glyphs);

   TeXBox& Add(const DArray<TeXBox>& tex_boxes);

   TeXBox& Add(const DArray<SPtr<TeXObject>>& tex_boxes);

   TeXBox& Add(TeXGlyph&& tex_glyph);

   TeXBox& Add(TeXBox&& tex_box);

   TeXBox& Add(DArray<TeXGlyph>&& tex_glyphs);

   TeXBox& Add(DArray<TeXBox>&& tex_boxes);

   TeXBox& SetPixelDensity(UInt value);

   TeXBox& SetAnchor(const SVectorR3& anchor);

   TeXBox& SetFontSize(const UChar font_size);

   TeXBox& SetName(const std::string& name) override;

   TeXBox& SetColour(const SVectorR4& rgba_colour) override;

   TeXBox& SetColour(const Colour& colour) override;

   TeXBox& SetItalic(bool is_italic) override;

   TeXBox& SetBold(bool is_bold) override;

 private:
   friend class Visualiser;

   void InitTeXBox(size_t id);

   void InitTeXObject(GlyphSheet::IndexT& index_offset) override;

   void ComputeDimensions(const GlyphSheet& glyph_sheet, UChar font_size, const SVectorR3& texbox_anchor, const SVectorR2& texbox_dimensions,
                          TeXSpacer& spacer) override;

   void OffsetAlongZ(Real& cumu_offset) override;

   fm::Path ImagePath() const;

   DArray<SPtr<TeXObject>> SubBoxes_;
   GlyphSheet              GlyphSheet_;
   SVectorR2               Dimensions_{}; // [width, height] in world-space coordinates.
   SVectorR3               Anchor_{};     // Bottom-left corner in world-space coordinates.
   UChar                   FontSize_{10}; // In units of LaTeX points (pt).

   /** Friendly unit tests */
   friend class ParseTeXTest_ParseTeXObject_Test;
};


}
