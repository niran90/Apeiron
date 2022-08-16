#include "../include/ParseTeX.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeX Parsing Functions
***************************************************************************************************************************************************************/
Glyph ParseTeXChar(const char c)
{
   Glyph glyph(c);
   if(c == OneOf(' ', '\t', '\n', '$')) glyph.DoNotRender();
   return glyph;
}

/***************************************************************************************************************************************************************
* TeX Parsing Helper Functions
***************************************************************************************************************************************************************/
bool isGlyphString(const std::string_view& tex_str)
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