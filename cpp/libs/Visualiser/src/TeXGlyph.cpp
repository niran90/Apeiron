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

#include "../include/TeXGlyph.h"
#include "../include/ModelFactory.h"
#include "../include/ParseTeX.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Glyph Public Interface
***************************************************************************************************************************************************************/
TeXGlyph::TeXGlyph(const char tex_char)
   : TeXGlyph(std::string(1, tex_char)) {}

TeXGlyph::TeXGlyph(const std::string& tex_str) { Set(tex_str); }

TeXGlyph&
TeXGlyph::Set(const char tex_char) { return Set(std::string(1, tex_char)); }

TeXGlyph&
TeXGlyph::Set(const std::string& tex_str)
{
   DEBUG_ASSERT(isGlyphString(tex_str), "The following string does not yet qualify as a glyph: ", tex_str)

   Text_ = tex_str;
   if(TeXSpacer::Required(tex_str))
   {
      AddSpacer_ = true;
      Text_ += TeXSpacer::Text();
   }

   return *this;
}

TeXGlyph&
TeXGlyph::SetColour(const SVectorR4& rgba_colour) { return SetColour(Colour{rgba_colour}); }

TeXGlyph&
TeXGlyph::SetColour(const Colour& colour)
{
   if(!Colour_) Colour_ = colour;
   return *this;
}

TeXGlyph&
TeXGlyph::SetItalic(const bool is_italic)
{
   if(!Italic_) Italic_ = is_italic;
   return *this;
}

TeXGlyph&
TeXGlyph::SetBold(const bool is_bold)
{
   if(!Bold_) Bold_ = is_bold;
   return *this;
}

/***************************************************************************************************************************************************************
* Glyph Private Interface
***************************************************************************************************************************************************************/
void
TeXGlyph::InitTeXObject(GlyphSheet::IndexT& index_offset)
{
   ASSERT(!Init_, "This glyph has already been initialised.")

   // Initialise glyph and assign glyph index if rendered.
   if(Rendered())
   {
      Init_  = true;
      Index_ = index_offset++;
   }
}

void
TeXGlyph::ComputeDimensions(const GlyphSheet& glyph_sheet, const UChar font_size, const SVectorR3& texbox_anchor, const SVectorR2& texbox_dims,
                            TeXSpacer& spacer)
{
   if(!Rendered()) return;

   ASSERT(Init_, "The glyph must be initialise before its dimensions are computed.")

   const auto& glyph_info = glyph_sheet.GlyphInfo(Index_);
   const auto  xy_scale   = GlyphSheet::FontSizeScale(font_size);
   const auto  y_scale    = 1.02; // Some glyphs slightly exceed the vertical bounds of the bounding box (e.g. m, e, 2), so need to vertically scale a bit.
   SVectorR2   glyph_dims;        // Dimensions of the glyph in the xy-plane in world-space coordinates.
   SVectorR2   glyph_anchor;      // Bottom-left corner of the glyph in world-space coordinates.

   // Compute the scaled glyph dimensions.
   glyph_dims.x() = xy_scale * static_cast<Real>(glyph_info.Width);
   glyph_dims.y() = xy_scale * y_scale * static_cast<Real>(glyph_info.Height + glyph_info.Depth);
   glyph_anchor   = xy_scale * SVectorR2{glyph_info.Position.x(), glyph_info.Position.y() - glyph_info.Depth }; // Anchor in TeXBox local coordinate system.

   // Ensure the glyph's vertical bounds haven't exceeded that of the tex-box after scaling.
   if(glyph_anchor.y() + glyph_dims.y() > texbox_dims.y()) Clip(glyph_dims.y(), Zero, texbox_dims.y() - glyph_anchor.y());

   // Set texture coordinates based on the glyph's dimensions w.r.t. the tex-box's dimensions.
   Mesh_ = ModelFactory::Rectangle(glyph_dims.x(), glyph_dims.y(), false).ModelMesh();
   auto set_tex_coor = [&](const size_t i, const SVectorR2& point)
      { Mesh_.Vertices_[i].TextureCoordinates = glm::vec2(point.x() / texbox_dims.x(), point.y() / texbox_dims.y()); };
   set_tex_coor(0, glyph_anchor);
   set_tex_coor(1, glyph_anchor + SVectorR2{glyph_dims.x(), Zero });
   set_tex_coor(2, glyph_anchor + glyph_dims);
   set_tex_coor(3, glyph_anchor + SVectorR2{Zero, glyph_dims.y() });

   // Apply cumulative spacer offsets scaled by the font size. Note: it is very important that this be added after texture coordinates are computed, since the
   // positive offsets (following characters requiring spacers) in the texture itself still remain.
   glyph_anchor.x() += font_size * spacer.Offset(glyph_anchor.x(), AddSpacer_);

   // Offset the glyph w.r.t. the parent tex-box's anchor. Also need to offset by half the glyph's dimensions as the glyph's rectangle is centred on the origin.
   glyph_anchor += Half * SVectorR2{glyph_dims.x(), glyph_dims.y() };
   glyph_anchor += SVectorR2{texbox_anchor.x(), texbox_anchor.y() };
   OffsetPosition(ToVector<3>(glyph_anchor));

   // Initialise underlying model.
   Model::Init();
}

void
TeXGlyph::OffsetAlongZ(Real& cumu_offset)
{
   this->OffsetPosition({0.0, 0.0, cumu_offset});
   cumu_offset += OffsetZ_;
}

}