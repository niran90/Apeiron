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
#include "TeXBox.h"
#include "TeXGlyph.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeX Parsing Functions
***************************************************************************************************************************************************************/
template<CharIterator It>
DArray<SPtr<TeXObject>>
ParseTeXString(It& current, const It last, const bool is_math_mode)
{
   DArray<SPtr<TeXObject>> tex_objects;
   tex_objects.reserve(std::distance(current, last));

   while(current < last)
   {
      const auto c = *current;
      if(c != '\\')
      {
         if(c != '$')           tex_objects.Append(ParseTeXObject(current, last, is_math_mode)); // TeX single/compound glyph
         else if(!is_math_mode) tex_objects.Append(ParseTeXMath(current, last));                 // TeX math string
         else EXIT("This code should not be reachable if already in math mode.")
      }
      else tex_objects.Append(ParseTeXCommand(current, last)); // TeX command
   }

   ASSERT(current == last, "Could not parse the TeX string correctly.")
   tex_objects.shrink_to_fit();
   return tex_objects;
}

template<CharIterator It>
SPtr<TeXObject>
ParseTeXObject(It& current, const It last, const bool is_math_mode)
{
   auto tex_glyph = ParseTeXChar(*(current++));

   // Math-mode and glyph is compound (i.e. contains trailing superscripts/subscripts).
   if(is_math_mode && tex_glyph->Rendered() && *current == OneOf('_', '^'))
   {
      // Accumulate trailing superscripts/subscripts.
      auto [sub_boxes, script_text] = ParseAllTeXScriptText(current, last);

      auto tex_box = std::make_shared<TeXBox>();
      tex_box->SetText(tex_glyph->Text() + script_text);
      tex_box->Add(tex_glyph);
      tex_box->Add(sub_boxes);

      return tex_box;
   }
   // Not math-mode OR no trailing superscripts/subscripts.
   else return tex_glyph;
}

template<CharIterator It>
SPtr<TeXObject>
ParseTeXCommand(It& current, const It last, const bool is_math_mode)
{
   ASSERT(*current == '\\', "A TeX command must begin with a backslash character.")

   const auto [cmd_prefix_end, has_args, render] = GetTeXCommandInfo(current, last);
   const auto command    = std::string(current, cmd_prefix_end);
   const auto trail_char = cmd_prefix_end != last ? *cmd_prefix_end : '\0';
   current = cmd_prefix_end;

   // Allocate a TeX object (TeXGlyph or TeXBox) depending on whether the command is to be rendered and whether it has arguments.
   SPtr<TeXObject> tex_object;
   if(render && has_args) tex_object = std::make_shared<TeXBox>();
   else
   {
      auto tex_glyph = std::make_shared<TeXGlyph>();
      if(!render) tex_glyph->DoNotRender();
      tex_object = std::move(tex_glyph);
   }

   tex_object->SetText(command);

   // Parse arguments (including trailing superscripts/subscripts) if the command has any.
   if(has_args)
   {
      // No trailing superscripts/subscripts, so parse {}-enclosed arguments.
      if(trail_char != OneOf('_', '^'))
      {
         ASSERT(trail_char == '{', "Expected at least one argument after the command ", command, ", but there is no trailing { brace.")

         // Get all enclosed arguments.
         const auto enclosures = GetFirstEnclosureChain(cmd_prefix_end, last, '{', '}', true);
         ASSERT(enclosures.size() == OneOf(1u, 2u), "Expected either one or two arguments for command ", command, ", but ", enclosures.size()," were found.")

         // Parse enclosures and (if rendered) add all enclosed objects as sub-objects of the parent tex-box.
         FOR_EACH(enclosure, enclosures)
         {
            // Parse sub-object/s in current enclosure.
            current = enclosure.first + 1;
            auto sub_boxes = ParseTeXString(current, enclosure.second - 1, is_math_mode);

            // Add opening brace, sub-glyph text, and closing brace.
            std::string encl_str;
            encl_str += '{';
            FOR_EACH(sub_box, sub_boxes) encl_str += sub_box->Text();
            encl_str += '}';

            // Update parent tex-box and add sub-objects, if rendered.
            tex_object->AddText(encl_str);
            if(render)
            {
               auto tex_box = std::dynamic_pointer_cast<TeXBox>(tex_object);
               ASSERT(tex_box, "Could not down-cast to a TeXBox.")
               tex_box->Add(sub_boxes);
            }
         }

         // Move one past the closing brace of the final enclosure.
         ++current;
      }
      // Has trailing superscripts/subscripts.
      else
      {
         ASSERT(render, "All LaTeX commands with trailing superscripts/subscripts must be rendered.")
         ASSERT(trail_char == OneOf('_', '^'),
                "Expected to parse a superscript/subscript but got the following trailing character: ", std::string{trail_char}, ".")

         // Downcast to tex-box to add sub-glyphs.
         auto tex_box = std::dynamic_pointer_cast<TeXBox>(tex_object);
         ASSERT(tex_box, "Could not down-cast to a TeXBox.")

         // Accumulate trailing superscripts/subscripts.
         auto [sub_boxes, script_text] = ParseAllTeXScriptText(current, last);
         tex_box->AddText(script_text);
         tex_box->Add(TeXGlyph(command));
         tex_box->Add(sub_boxes);
      }
   }

   return tex_object;
}

template<CharIterator It>
DArray<SPtr<TeXObject>>
ParseTeXMath(It& current, const It last)
{
   ASSERT(*current == '$', "In math mode, the first character must point to the $ symbol.")

   DArray<SPtr<TeXObject>> tex_objects;
   tex_objects.reserve(std::distance(current, last));

   const auto brace_glyph = ParseTeXChar('$');
   const bool is_double_braced = *(current + 1) == '$';

   // Add opening math-mode brace/s.
   tex_objects.Append(brace_glyph);
   if(is_double_braced) tex_objects.Append(brace_glyph);

   // Find the bounding iterators of the enclosing math-mode string.
   if(is_double_braced) ++current;
   const auto [bounds, found] = GetFirstEnclosure(current, last, '$', '$');
   ASSERT(found, "Could not find an enclosed substring in math mode.")

   // Add glyphs from the enclosing string.
   current = bounds.first;
   tex_objects.Append(ParseTeXString(current, bounds.second, true));
   ASSERT(current == bounds.second, "Could not correctly parse the enclosed substring in math mode.")

   // Add closing math-mode brace/s.
   tex_objects.Append(brace_glyph);
   ++current;
   if(is_double_braced)
   {
      ASSERT(*(bounds.second + 1) == '$', "Expected a second $ symbol, but did not find one.")
      tex_objects.Append(brace_glyph);
      ++current;
   }

   tex_objects.shrink_to_fit();
   return tex_objects;
}

template<CharIterator It>
std::pair<DArray<SPtr<TeXObject>>, std::string>
ParseAllTeXScriptText(It& current, const It last)
{
   std::string script_str;
   DArray<SPtr<TeXObject>> tex_objects;
   tex_objects.reserve(2);

   FOR(i, 2)
      if(*current == OneOf('_', '^'))
      {
         const bool is_braced = *(current + 1) == '{';

         // Add underscore/carat and opening brace for sub-glyph insertion.
         script_str += *current;
         if(is_braced) script_str += '{';

         // Add subglyphs.
         auto script_glyphs = ParseTeXScriptText(current, last);
         FOR_EACH(glyph, script_glyphs) script_str += glyph->Text();

         // Add closing brace and append subglyphs.
         if(is_braced) script_str += '}';
         tex_objects.Append(std::move(script_glyphs));
      }

   tex_objects.shrink_to_fit();
   return { tex_objects, script_str };
}

/***************************************************************************************************************************************************************
* TeX Parsing Helper Functions
***************************************************************************************************************************************************************/
template<CharIterator It>
bool
isTeXCommandEnd(const It current, const It last)
{
   if(current == last) return true;
   const auto c = *current;
   return c == OneOf(' ', '.', ',', ';', ':', '(', '[', '\t', '\n', '\\') || std::isdigit(c) ||
          (!std::isalpha(*(current - 1)) && std::isalpha(c)); // Note: assumes base string has at least one char before current iterator.
}

template<CharIterator It>
std::pair<bool, It>
isTeXCharCommand(const It first, const It last)
{
   ASSERT(*first == '\\', "A TeX command must begin with a backslash character.")

   const bool is_char_cmd = std::distance(first, last) > 1 && !std::isdigit(*(first + 1)) &&
                           (std::distance(first, last) > 2 ? isTeXCommandEnd(first + 2, last) : true);
   const auto end_iter = is_char_cmd ? GetTeXCommandPrefixEnd(first, last, true) : last;

   return { is_char_cmd, end_iter };
}

template<CharIterator It>
std::pair<bool, It>
isTeXWordCommand(const It first, const It last)
{
   ASSERT(*first == '\\', "A TeX command must begin with a backslash character.")

   const bool is_word_cmd = std::distance(first, last) > 2 &&
                            std::all_of(first + 1, first + 3, [](char c){ return std::isalpha(c); }); // At least two letters in the command name
   const auto end_iter = is_word_cmd ? GetTeXCommandPrefixEnd(first, last, false) : last;

   return { is_word_cmd, end_iter };
}

template<CharIterator It>
std::tuple<It, bool, bool>
GetTeXCommandInfo(const It current, const It last)
{
   ASSERT(*current == '\\', "A TeX command must begin with a backslash character.")
   ASSERT(std::distance(current, last) > 1, "A TeX command must comprise at least two characters.")

   bool char_cmd{}, word_cmd{};
   It cmd_prefix_end;
   std::tie(char_cmd, cmd_prefix_end)               = isTeXCharCommand(current, last);
   if(!char_cmd) std::tie(word_cmd, cmd_prefix_end) = isTeXWordCommand(current, last);
   ASSERT(char_cmd != word_cmd, "Expected either a command character or a command word, but got ", std::string(current, last), ".")

   const auto cmd_prefix = std::string(current + 1, cmd_prefix_end); // Command name without the backslash

   if(char_cmd) return { cmd_prefix_end, false, cmd_prefix != OneOf("!", ",", ":", ">", ";", " ") }; // Single characters with no trailing args.
                                                                                                     // Render if not a spacing character.
   else if(isTeXCommandEnd(cmd_prefix_end, last)) return { cmd_prefix_end, false, cmd_prefix != OneOf("quad", "qquad", "linebreak") }; // Word command without trailing arguments.
                                                                                                                                       // Render if not a spacing word.
   else return { cmd_prefix_end, true, cmd_prefix != OneOf( "dec",  "decc",  "deccc",  "inc",  "incc",  "inccc",
                                                           "vdec", "vdecc", "vdeccc", "vinc", "vincc", "vinccc") }; // Word command with trailing arguments.
                                                                                                                    // Render if not a spacing word.
}

template<CharIterator It>
DArray<SPtr<TeXObject>>
ParseTeXScriptText(It& current, const It last)
{
   ASSERT(*current == OneOf('_', '^'), "TeX superscript/subscript text must begin with either an underscore or a carat.")

   DArray<SPtr<TeXObject>> tex_objects;
   tex_objects.reserve(std::distance(current, last));

   const auto next = current + 1;
   ASSERT(next != last, "Incomplete TeX string.")
   const auto next_char = *next;

   // Enclosed string
   if(next_char == '{')
   {
      const auto [bounds, found] = GetFirstEnclosure(current, last, '{', '}');
      ASSERT(found, "Could not find an enclosed substring in the subscript/superscript.")

      // Add glyphs from the enclosed string.
      current = bounds.first;
      tex_objects.Append(ParseTeXString(current, bounds.second, true));
      ASSERT(current == bounds.second, "Could not correctly parse the enclosed substring in the subscript/superscript.")
      ++current; // Advance one past the closing brace
   }
   // Single character
   else
   {
      tex_objects.Append(std::make_shared<TeXGlyph>(next_char));
      current += 2;
   }

   tex_objects.shrink_to_fit();
   return tex_objects;
}

template<CharIterator It>
It GetTeXCommandPrefixEnd(const It first, const It last, const bool is_char_cmd)
{
   DEBUG_ASSERT(*first == '\\', "A TeX command must begin with a backslash character.")
   return is_char_cmd ? first + 2 : std::find_if_not(first + 1, last, [](auto c){ return std::isalpha(c); });
}

}