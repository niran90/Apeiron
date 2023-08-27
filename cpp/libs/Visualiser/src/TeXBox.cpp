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

#include "../include/ParseTeX.h"
#include "../include/TeXBox.h"
#include "../include/Texture.h"
#include "FileManager/include/File.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeXBox Public Interface
***************************************************************************************************************************************************************/
TeXBox::TeXBox(const std::string& str) { Add(str); }

TeXBox::TeXBox(const DArray<TeXGlyph>& tex_glyphs) { Add(tex_glyphs); }

TeXBox::TeXBox(const DArray<TeXBox>& tex_boxes) { Add(tex_boxes); }

TeXBox::TeXBox(DArray<TeXGlyph>&& tex_glyphs) { Add(std::move(tex_glyphs)); }

TeXBox::TeXBox(DArray<TeXBox>&& tex_boxes) { Add(std::move(tex_boxes)); }

TeXBox&
TeXBox::Add(const std::string& str) { return Add(ParseTeXString(str)); }

TeXBox&
TeXBox::Add(const TeXGlyph& tex_glyph)
{
   auto glyph = std::make_shared<TeXGlyph>(tex_glyph);
   SubBoxes_.push_back(glyph);
   SubModels_.push_back(glyph);
   return *this;
}

TeXBox&
TeXBox::Add(const TeXBox& tex_box)
{
   auto box = std::make_shared<TeXBox>(tex_box);
   SubBoxes_.push_back(box);
   SubModels_.push_back(box);
   return *this;
}

TeXBox&
TeXBox::Add(const SPtr<TeXObject>& tex_object)
{
   // Add as a sub-box.
   SubBoxes_.push_back(tex_object);

   // Down-cast to either a TeXBox or a TeXGlyph.
   auto tex_box   = std::dynamic_pointer_cast<TeXBox>(tex_object);
   auto tex_glyph = std::dynamic_pointer_cast<TeXGlyph>(tex_object);
   ASSERT(tex_box || tex_glyph, "Failed to down-cast to either a TeXBox or a TeXGlyph.")

   // Add down-casted object as a sub-model.
   if(tex_box) SubModels_.push_back(tex_box);
   else        SubModels_.push_back(tex_glyph);
   return *this;
}

TeXBox&
TeXBox::Add(const DArray<TeXGlyph>& tex_glyphs)
{
   FOR_EACH(tex_glyph, tex_glyphs) Add(tex_glyph);
   return *this;
}

TeXBox&
TeXBox::Add(const DArray<TeXBox>& tex_boxes)
{
   FOR_EACH(tex_box, tex_boxes) Add(tex_box);
   return *this;
}

TeXBox&
TeXBox::Add(const DArray<SPtr<TeXObject>>& tex_boxes)
{
   FOR_EACH(tex_box, tex_boxes) Add(tex_box);
   return *this;
}

TeXBox&
TeXBox::Add(TeXGlyph&& tex_glyph)
{
   auto glyph = std::make_shared<TeXGlyph>(std::move(tex_glyph));
   SubBoxes_.push_back(glyph);
   SubModels_.push_back(glyph);
   return *this;
}

TeXBox&
TeXBox::Add(TeXBox&& tex_box)
{
   auto box = std::make_shared<TeXBox>(std::move(tex_box));
   SubBoxes_.push_back(box);
   SubModels_.push_back(box);
   return *this;
}

TeXBox&
TeXBox::Add(DArray<TeXGlyph>&& tex_glyphs)
{
   FOR_EACH(tex_glyph, tex_glyphs) Add(std::move(tex_glyph));
   return *this;
}

TeXBox&
TeXBox::Add(DArray<TeXBox>&& tex_boxes)
{
   FOR_EACH(tex_box, tex_boxes) Add(std::move(tex_box));
   return *this;
}

TeXBox&
TeXBox::SetPixelDensity(const UInt value)
{
   GlyphSheet_.SetPixelDensity(value);
   return *this;
}

TeXBox&
TeXBox::SetAnchor(const SVectorR3& anchor)
{
   Anchor_ = anchor;
   return *this;
}

TeXBox&
TeXBox::SetFontSize(const UChar font_size)
{
   FontSize_ = font_size;
   return *this;
}

TeXBox&
TeXBox::SetName(const std::string& name)
{
   ModelGroup::SetName(name);
   return *this;
}

TeXBox&
TeXBox::SetColour(const SVectorR4& rgba_colour) { return SetColour(Colour{rgba_colour}); }

TeXBox&
TeXBox::SetColour(const Colour& colour)
{
   FOR_EACH(str, SubBoxes_) str->SetColour(colour);
   return *this;
}

TeXBox&
TeXBox::SetItalic(const bool is_italic)
{
   FOR_EACH(str, SubBoxes_) str->SetItalic(is_italic);
   return *this;
}

TeXBox&
TeXBox::SetBold(const bool is_bold)
{
   FOR_EACH(str, SubBoxes_) str->SetBold(is_bold);
   return *this;
}

/***************************************************************************************************************************************************************
* TeXBox Private Interface
***************************************************************************************************************************************************************/
void
TeXBox::InitTeXBox(const size_t id)
{
   // Propagate name (if already set) to sub-glyphs.
   if(!Name_.empty()) ModelGroup::SetName(Name_);

   // Initialise sub-boxes, assign glyph indices, and accumulate text.
   GlyphSheet::IndexT glyph_index = 0;
   InitTeXObject(glyph_index);

   // Initialise glyph sheet.
   GlyphSheet_.Init(id, Text_);
   ASSERT(GlyphSheet_.Width() && GlyphSheet_.Height(), "The dimensions of the glyph sheet must be computed before those of the TeXBox.")

   // Compute the world-space dimensions from the glyph-sheet dimensions, which are currently expressed in scaled points.
   Dimensions_ = { GlyphSheet_.Width(), GlyphSheet_.Height() };
   Dimensions_ *= GlyphSheet::FontSizeScale(FontSize_); // Scaled points -> points -> world-space * font-size.

   // Compute the sub-glyph dimensions and their texture coordinates.
   TeXSpacer spacer;
   ComputeDimensions(GlyphSheet_, FontSize_, Anchor_, Dimensions_, spacer);

   // Offset along the z-axis to avoid z-fighting.
   Real offset{};
   OffsetAlongZ(offset);
}

void
TeXBox::InitTeXObject(GlyphSheet::IndexT& index_offset)
{
   // Initialise sub-boxes and add contributions to the TeX-box string, if it has not already been set.
   const auto set_text = Text_.empty();
   FOR_EACH(sub_box, SubBoxes_)
   {
      sub_box->InitTeXObject(index_offset);
      if(set_text) Text_ += sub_box->Text();
   }
}

void
TeXBox::ComputeDimensions(const GlyphSheet& glyph_sheet, const UChar font_size, const SVectorR3& texbox_anchor, const SVectorR2& texbox_dimensions,
                          TeXSpacer& spacer)
{
   FOR_EACH(sub_box, SubBoxes_)
      sub_box->ComputeDimensions(glyph_sheet, font_size, texbox_anchor, texbox_dimensions, spacer);
}

fm::Path
TeXBox::ImagePath() const
{
   const auto compile_dir = GlyphSheet_.CompileDirectory();
   ASSERT(!compile_dir.empty(), "The compile directory has not yet been set for this TeXBox.")
   return compile_dir / LaTeXTemplate().filename().replace_extension(".png");
}

void
TeXBox::OffsetAlongZ(Real& cumu_offset)
{
   FOR_EACH(sub_box, SubBoxes_)
      sub_box->OffsetAlongZ(cumu_offset);
}

}