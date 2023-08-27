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

#include "../include/GlyphSheet.h"
#include "FileManager/include/File.h"

namespace aprn::vis {

void
GlyphSheet::Init(size_t id, const std::string& text)
{
   ASSERT(!text.empty(), "The TeXBox being compiled has no text.")

   // Set the compile-directory for this glyph sheet.
   CompileDirectory_ = LaTeXDirectory() / ("texbox" + ToString(id));

   // Compile LaTeX and generate image.
   CompileLaTeXSource(text);
   CreateGlyphSheetImage();

   // Read glyph data and size the glyph sheet.
   ReadGlyphBoxPositions();
   ReadGlyphBoxAttributes();
   ComputeDimensions();
}

void
GlyphSheet::CompileLaTeXSource(const std::string& text)
{
   const auto& comp_dir = CompileDirectory_;

   // Initialise the LaTeX compile directory for this TeX-box and copy over the LaTeX template.
   fm::CreateDirectory(comp_dir);
   fm::ClearDirectory(comp_dir);
   fm::CopyFile(LaTeXTemplate() , comp_dir);
   fm::CopyFile(LuaTeXTemplate(), comp_dir);

   // Transfer TeX-box text to the .tex file.
   TeXFile_ = comp_dir / LaTeXTemplate().filename();
   fm::File file(TeXFile_, fm::Mode::Append);
   file.Write("\n", text, "\n\\end{document}");
   file.Close();

   // Compile LaTeX source code.
   fm::CompileTeXFile("lualatex", TeXFile_);
}

void
GlyphSheet::CreateGlyphSheetImage() { fm::ConvertPDFtoPNG(TeXFile_.replace_extension(".pdf"), PixelDensity_); }

void
GlyphSheet::ReadGlyphBoxPositions()
{
   fm::File file(CompileDirectory_ / "positions.txt", fm::Mode::Read);
   IndexT glyph_index{};

   while(!file.isEnd())
   {
      Boxes_.emplace_back();
      auto& x = Boxes_[glyph_index].Position.x();
      auto& y = Boxes_[glyph_index].Position.y();

      file.Read(x, y);
      if(x < 0)
      {
         DEBUG_ASSERT(x == -1, "Read in an x-coordinate which is not -1: ", x)
         x = 0;
      }
      ++glyph_index;
   }
   file.Close();
}

void
GlyphSheet::ReadGlyphBoxAttributes()
{
   // Read glyph box attributes. Note: need to read with a wide file, as the glyph characters must be read in as wchar_t.
   fm::WFile wfile;
   wfile.Open(CompileDirectory_ / "attributes.txt", fm::Mode::Read);
   IndexT glyph_index{};

   while(wfile.isValid())
   {
      DEBUG_ASSERT(glyph_index < Boxes_.size(), "The glyph index ", glyph_index, " exceeds the total number of glyphs.")

      auto& c = Boxes_[glyph_index].Char;
      auto& w = Boxes_[glyph_index].Width;
      auto& h = Boxes_[glyph_index].Height;
      auto& d = Boxes_[glyph_index].Depth;

      wfile.Read(c, w, h, d);
      ++glyph_index;
   }

   wfile.Close();
   ASSERT(glyph_index == Boxes_.size(), "The number of glyph attributes does not match the number of positions read in.")
}

void
GlyphSheet::ComputeDimensions()
{
   using int_T  = decltype(Width_);
   using coor_T = SVector2<int_T>;

   // Calculate glyph sheet dimensions (min/max bounds of all glyph boxes).
   coor_T min_pos(MaxInt<int_T>);
   coor_T max_pos(MinInt<int_T>);

   FOR_EACH_CONST(glyph, Boxes_)
   {
      const coor_T bott_left = { glyph.Position.x()              , glyph.Position.y() - glyph.Depth  };
      const coor_T top_right = { glyph.Position.x() + glyph.Width, glyph.Position.y() + glyph.Height };

      FOR(i, 2)
      {
         min_pos[i] = Min(min_pos[i], bott_left[i]);
         max_pos[i] = Max(max_pos[i], top_right[i]);
      }
   }
   Width_  = max_pos.x() - min_pos.x();
   Height_ = max_pos.y() - min_pos.y();
}

const GlyphBox&
GlyphSheet::GlyphInfo(const GlyphSheet::IndexT glyph_index) const
{
   ASSERT(glyph_index < Boxes_.size(), "The glyph index ", glyph_index, " is out of bounds.")
   return Boxes_[glyph_index];
}

}
