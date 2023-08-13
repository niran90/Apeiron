/***************************************************************************************************************************************************************
* GPL-3.0 License
* Copyright (C) 2022 Niran A. Ilangakoon
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program.
* If not, see <https://www.gnu.org/licenses/>.
***************************************************************************************************************************************************************/

#pragma once

#include "../../../include/Global.h"
#include "DataContainer/include/Array.h"
#include "FileManager/include/FileSystem.h"
#include "Visualiser/include/Texture.h"

#include <unordered_map>

namespace aprn::vis {

namespace fm = flmgr;

/***************************************************************************************************************************************************************
* Glyph Box Struct
***************************************************************************************************************************************************************/
struct GlyphBox
{
   wchar_t         Char{};
   Int64           Width{};
   Int64           Height{};
   Int64           Depth{};
   SVector2<Int64> Position{};
};

/***************************************************************************************************************************************************************
* Glyph Sheet Class
***************************************************************************************************************************************************************/
class GlyphSheet
{
 public:
   typedef UInt16 IndexT;

   void Init(size_t id, const std::string& text);

   void CompileLaTeXSource(const std::string& text);

   void CreateGlyphSheetImage();

   void ReadGlyphBoxPositions();

   void ReadGlyphBoxAttributes();

   void ComputeDimensions();

   const GlyphBox& GlyphInfo(const IndexT glyph_index) const;

   inline auto Width() const { return Width_; }

   inline auto Height() const { return Height_; }

   inline auto CompileDirectory() const { return CompileDirectory_; }

   inline void SetPixelDensity(const UInt density) { PixelDensity_ = density; }

   constexpr static Real PointSize() { return PointSize_; }

   constexpr static Real ScaledPointsToPointsFactor() { return One / 65536.0; }

   constexpr static Real FontSizeScale(const UChar font_size) { return ScaledPointsToPointsFactor() * PointSize_ * static_cast<Real>(font_size); }

 private:
   DArray<GlyphBox>      Boxes_;
   fm::Path              CompileDirectory_;
   fm::Path              TeXFile_;
   Int64                 Width_{};            // Width in units of LaTeX scaled points (sp).
   Int64                 Height_{};           // Height in units of LaTeX scaled points (sp).
   UInt16                PixelDensity_{3000}; // Measured in DPI.
   constexpr static Real PointSize_{0.005};   // Height of 1pt (in LaTeX) in world-space.
   constexpr static bool ShowOutput_{false};  // Show output of LaTeX compilation and pdf->png image conversion.
};

/***************************************************************************************************************************************************************
* Stand-alone Functions
***************************************************************************************************************************************************************/
inline fm::Path LaTeXDirectory() { return "./libs/Visualiser/data/latex"; }

inline fm::Path LaTeXTemplate() { return "./libs/Visualiser/resources/latex/texbox.tex"; }

inline fm::Path LuaTeXTemplate() { return "./libs/Visualiser/resources/latex/write_boxes.lua"; }

inline void InitTeXDirectory() { fm::CreateDirectory(LaTeXDirectory(), true); }

}