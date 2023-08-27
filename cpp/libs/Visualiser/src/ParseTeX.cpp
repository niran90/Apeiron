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

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeX Parsing Functions
***************************************************************************************************************************************************************/
SPtr<TeXGlyph>
ParseTeXChar(const char c)
{
   auto glyph = std::make_shared<TeXGlyph>(c);
   if(c == OneOf(' ', '\t', '\n', '$')) glyph->DoNotRender();
   return glyph;
}

DArray<SPtr<TeXObject>>
ParseTeXString(const std::string& str)
{
   auto current = str.begin();
   return ParseTeXString(current, str.end());
}

/***************************************************************************************************************************************************************
* TeX Parsing Helper Functions
***************************************************************************************************************************************************************/
bool
isGlyphString(const std::string_view& tex_str)
{
   if(tex_str.empty()) return false;

   const bool is_cmd = tex_str.front() == '\\';

   if(tex_str.length() == (!is_cmd ? 1 : 2)) return true; // Single character glyphs/TeX commands
   else if(is_cmd) // All TeX word commands
   {
      std::string bare_str(tex_str);
      const auto [is_word, _] = isTeXWordCommand(bare_str.begin(), bare_str.end());

      if(is_word)
      {
         // First remove all chained enclosures, if this command has any.
         const auto enclosure_chain = GetAllEnclosures(bare_str, '{', '}', true);
         FOR_EACH_CONST(enclosure, enclosure_chain) bare_str = Remove(enclosure, bare_str);

         // Must only compute end iterator here as it may have been invalidated after enclosure removal.
         const auto   end_iter       = GetTeXCommandPrefixEnd(bare_str.begin(), bare_str.end(), false);
         const size_t cmd_prefix_len = std::distance(bare_str.begin(), end_iter);
         ASSERT(cmd_prefix_len <= bare_str.length(), "Could not correctly identify the end of the command prefix.")

         return cmd_prefix_len == bare_str.length() || std::all_of(end_iter, bare_str.end(), [](char c){ return c == OneOf('_', '^'); }); // e.g. \sum_{i=0}^N
      }
   }
   return false;
}

}