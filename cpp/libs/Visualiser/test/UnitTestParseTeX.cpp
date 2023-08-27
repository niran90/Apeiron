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

#include <gtest/gtest.h>
#include "../include/ParseTeX.h"
#include "../include/Window.h"

#ifdef DEBUG_MODE

namespace aprn::vis {

class ParseTeXTest : public testing::Test
{
 public:
   ParseTeXTest() {}

   void
   SetUp() override { TeXSpacer::Disable(); }
};

/***************************************************************************************************************************************************************
* Test TeX Parsing Functions
***************************************************************************************************************************************************************/
TEST_F(ParseTeXTest, ParseTeXChar)
{
   for(const auto c: {'.', ',', ';', ':', '`', '"', '(', '[', '{', '1', 'a', '?', '!'})
   {
      const auto glyph = ParseTeXChar(c);
      EXPECT_EQ(glyph->Text(), std::string(1, c));
      EXPECT_TRUE(glyph->Rendered());
   }
   for(const auto c: {' ', '\t', '\n', '$'})
   {
      const auto glyph = ParseTeXChar(c);
      EXPECT_EQ(glyph->Text(), std::string(1, c));
      EXPECT_FALSE(glyph->Rendered());
   }
   EXPECT_DEATH(ParseTeXChar('\\'), "");
}

TEST_F(ParseTeXTest, ParseTeXString)
{
   for(const std::string tex_str: {R"(The final velocity at time $t$ is computed as $v = u + at$.)",
                                   R"(Einstein showed that $e = mc^2$.)",
                                   R"(Euler's identity can be expressed as $e^{i\pi} = 0$.)",
                                   R"(Useful binomial identity: $$(a - b)^2 = (a + b)(a - b)$$.)",
                                   R"(An n\textsuperscript{th}-order polynomial can be expressed as $$ P_n(x) = \sum_{i = 0}^n c_ix^i. $$)",
                                   R"(The Divergence Theorem for volume $V$ bounded by surface $A$: $\int_V \nabla\cdot\bm{u} dV = \oint_A \bm{u}\cdot\bm{n} dA$.)"})
   {
      auto iter = tex_str.begin();
      const auto glyphs = ParseTeXString(iter, tex_str.end());

      std::string glyph_str;
      FOR_EACH(glyph, glyphs) glyph_str.append(glyph->Text());
      EXPECT_EQ(glyph_str, tex_str);
      EXPECT_EQ(iter, tex_str.end());
   }
}

TEST_F(ParseTeXTest, ParseTeXObject)
{
   // Not math mode
   for(const std::string g : {"a", "A", "(", "{", "[", ",", ":"})
   {
      auto iter = g.begin();
      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(ParseTeXObject(iter, g.end(), false));

      EXPECT_TRUE(glyph);
      EXPECT_EQ(glyph->Text(), g);
      EXPECT_EQ(iter, g.end());
      EXPECT_TRUE(glyph->Rendered());
   }
   for(const std::string g : {" ", "\t", "\n", "$"})
   {
      auto iter = g.begin();
      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(ParseTeXObject(iter, g.end(), false));

      EXPECT_TRUE(glyph);
      EXPECT_EQ(glyph->Text(), g);
      EXPECT_EQ(iter, g.end());
      EXPECT_FALSE(glyph->Rendered());
   }

   // Math mode - compound glyphs.
   for(const std::string g : {"x^2 ", "y_i^{3},", "u_t.", "v^{2}_i:"})
   {
      auto iter = g.begin();
      const auto tex_box = std::dynamic_pointer_cast<TeXBox>(ParseTeXObject(iter, g.end(), true));
      EXPECT_TRUE(tex_box);
      const auto script_count = std::count_if(g.begin(), g.end(), [](auto c){ return c == OneOf('_', '^'); });

      std::string str_check = g;
      std::string str       = tex_box->Text();
      str_check = Remove("{", str_check);
      str_check = Remove("}", str_check);
      str       = Remove("{", str);
      str       = Remove("}", str);

      EXPECT_EQ(str_check.substr(0, str.length()), str);
      EXPECT_EQ(tex_box->SubBoxes_.size(), script_count + static_cast<bool>(script_count));
      const auto iter_check = std::find_if(g.begin(), g.end(), [](auto c){ return c == OneOf(' ', ',', '.', ':'); });
      EXPECT_EQ(iter, iter_check);
   }
}

TEST_F(ParseTeXTest, ParseTeXCommand)
{
   // Characters
   for(const std::string cmd: {"\\p", "\\%,", "\\$: ", "\\^, ", "\\#:"})
   {
      auto iter = cmd.begin();
      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(ParseTeXCommand(iter, cmd.end(), false));
      EXPECT_TRUE(glyph);

      EXPECT_EQ(glyph->Text(), cmd.substr(0, 2));
      EXPECT_EQ(iter, cmd.begin() + 2);
      EXPECT_TRUE(glyph->Rendered());
   }

   // Characters (spacing)
   for(const std::string cmd: {"\\ ", "\\!", "\\, ", "\\;end", "\\ end"})
   {
      auto iter = cmd.begin();
      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(ParseTeXCommand(iter, cmd.end(), false));
      EXPECT_TRUE(glyph);

      EXPECT_EQ(glyph->Text(), cmd.substr(0, 2));
      EXPECT_EQ(iter, cmd.begin() + 2);
      EXPECT_FALSE(glyph->Rendered());
   }

   auto find_end = [](const std::string& cmd)
      { return std::find_if(cmd.begin(), cmd.end(), [](auto c){ return c == OneOf(' ', ',', '.', ';', ':', '('); }); };

   // Non-math mode words
   for(const std::string cmd: {"\\linebreak", "\\textit{italic}", "\\textbf{bold} font", "\\cite{Newton}, "})
   {
      auto iter = cmd.begin();
      const auto object = ParseTeXCommand(iter, cmd.end(), false);
      const auto iter_check = find_end(cmd);
      EXPECT_EQ(object->Text(), std::string(cmd.begin(), iter_check));
      EXPECT_EQ(iter, iter_check);
   }

   // Non-math mode words (spacing)
   for(const std::string cmd: {"\\quad", "\\qquad", "\\quad,", "\\qquad ", "\\quad 5x", "\\qquad x = 1.",
                               "\\dec{1}" , "\\decc{1} " , "\\deccc{1}, ", "\\inc{1}: "   , "\\incc{1};" , "\\inccc{1}(x)",
                               "\\vdec{1}", "\\vdecc{1}.", "\\vdeccc{1} ", "\\vinc{1} = 5", "\\vincc{1} ", "\\vinccc{1}"})
   {
      auto iter = cmd.begin();
      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(ParseTeXCommand(iter, cmd.end(), false));
      EXPECT_TRUE(glyph);

      const auto iter_check = find_end(cmd);
      EXPECT_EQ(glyph->Text(), std::string(cmd.begin(), iter_check));
      EXPECT_EQ(iter, iter_check);
      EXPECT_FALSE(glyph->Rendered());
   }

   // Math mode
   for(const std::string cmd: {"\\pi.", "\\eta,", "\\tau^{2}:", "\\frac{x^{2} - 1}{x + 1} =",
                               "\\sum_{i = 0}^{N} x_i", "\\dotp{\\bm{u}}{\\bm{v}},", "\\frac{\\bm{u} \\cross \\bm{v}}{x + 1},"})
   {
      auto iter = cmd.begin();
      const auto [cmd_prefix_end, has_args, render] = GetTeXCommandInfo(iter, cmd.end());

      const auto iter_check_reverse = std::find_if(cmd.rbegin(), cmd.rend(), [](auto c){ return c == OneOf(' ', ',', '.', ':'); });
      const auto iter_check = (iter_check_reverse + 1).base();

      const auto object = ParseTeXCommand(iter, cmd.end(), true);
      EXPECT_EQ(iter, iter_check);
      EXPECT_EQ(object->Text(), std::string(cmd.begin(), iter_check));

      if(!(render && has_args)) EXPECT_TRUE(std::dynamic_pointer_cast<TeXGlyph>(object)->Rendered());
   }

   // Death tests
   for(const std::string cmd: {" ", "2", "pi", "p_1", "frac{1}{2}"})
   {
      auto iter = cmd.begin();
      EXPECT_DEATH(ParseTeXCommand(iter, cmd.end()), "");
   }
}

TEST_F(ParseTeXTest, ParseTeXMath)
{
   // Math mode expressions.
   for(const std::string tex_str: {R"($v = u + at.$ )",
                                   R"($e = mc^2, $ )",
                                   R"($e^{i\pi} = 0; $ )",
                                   R"($$(a - b)^2 = (a + b)(a - b)$$.)",
                                   R"($$ P_5(x) = c_0 + c_1x + c_2x^2 + c_3x^3 + c_4x^4 + c_5x^5 $$;)",
                                   R"($$ \int_V \nabla\cdot\bm{u}\ dV = \oint_A \bm{u}\cdot\bm{n} dA.$$,)"})
   {
      auto iter = tex_str.begin();
      const auto glyphs = ParseTeXMath(iter, tex_str.end());

      std::string glyph_str;
      FOR_EACH(glyph, glyphs) glyph_str.append(glyph->Text());
      EXPECT_EQ(glyph_str, tex_str.substr(0, glyph_str.length()));
      auto iter_check = std::find_if(tex_str.rbegin(), tex_str.rend(), [](auto c){ return c == OneOf(' ', ',', '.', ';'); });
      EXPECT_EQ(iter, (iter_check + 1).base());
   }

   // Death tests.
   for(const std::string tex_str: {R"(v = u + at.)",
                                   R"($e = mc^2)",
                                   R"($$(a - b)^2 = (a + b)(a - b)$)"})
   {
      auto iter = tex_str.begin();
      EXPECT_DEATH(ParseTeXMath(iter, tex_str.end()), "");
   }
}

TEST_F(ParseTeXTest, ParseAllTeXScriptText)
{
   for(const std::string str: {"x^2_{i}", "\\tau_0^3 ", "\\sum_{i}^N x", "\\int_{0}^T t^2 dt"})
   {
      auto iter_start = std::find_if(str.begin(), str.end(), [](auto c){return c == OneOf('_', '^');});
      auto iter = iter_start;
      auto [glyphs, script_text] = ParseAllTeXScriptText(iter, str.end());

      EXPECT_EQ(glyphs.size(), 2);
      EXPECT_EQ(std::string(iter_start, str.end()).substr(0, script_text.length()), script_text);

      const auto glyph0 = std::dynamic_pointer_cast<TeXGlyph>(glyphs[0]);
      const auto glyph1 = std::dynamic_pointer_cast<TeXGlyph>(glyphs[1]);
      EXPECT_TRUE(glyph0);
      EXPECT_TRUE(glyph1);
      EXPECT_TRUE(glyph0->Rendered());
      EXPECT_TRUE(glyph1->Rendered());

      // First script text
      if(*(iter_start + 1) != '{') EXPECT_EQ(glyph0->Text(), std::string(1, *(iter_start + 1))); // Not braced
      else EXPECT_EQ(glyph0->Text(), std::string(1, *(iter_start + 2))); // Braced

      // Second script text
      iter_start = std::find_if(iter_start + 1, str.end(), [](auto c){return c == OneOf('_', '^');});
      if(*(iter_start + 1) != '{') // Not braced
      {
         EXPECT_EQ(glyph1->Text(), std::string(1, *(iter_start + 1)));
         EXPECT_EQ(iter, iter_start + 2);
      }
      else // Braced
      {
         EXPECT_EQ(glyph1->Text(), std::string(1, *(iter_start + 2)));
         EXPECT_EQ(iter, iter_start + 4);
      }
   }
}

/***************************************************************************************************************************************************************
* Test TeX Parsing Helper Functions
***************************************************************************************************************************************************************/
TEST_F(ParseTeXTest, isGlyphString)
{
   // Test single characters and TeX commands (ordinary strings and raw literal strings)
   for(const auto c: {"a", "A", " ", "\t", "\n", "$", "(", "{", "[", ",", ":", "\\(", "\\{", "\\[", "\\,"}) EXPECT_TRUE(isGlyphString(c));
   for(const auto c: {"\\;", "\\tau", "\\frac{1}{2}", "\\frac{x^2 + 1}{x^2 - 1}", "\\sum_{i = 0}^{\\infty}"}) EXPECT_TRUE(isGlyphString(c));
   for(const auto c: {R"(\;)", R"(\tau)", R"(\frac{1}{2})", R"(\sum_{i = 0}^{\infty})", R"( )"}) EXPECT_TRUE(isGlyphString(c));

   for(const auto c: {"ab", "AB", "  ", "\tt", "\nn", "$$", "{}", "\\()", "\\[1]", "\\,a"}) EXPECT_FALSE(isGlyphString(c));
   for(const auto c: {"\\;1", "\\tau ", "\\frac{1}{2},", "\\sum_{i = 0}^{\\infty}\\pi"}) EXPECT_FALSE(isGlyphString(c));
   for(const auto c: {R"(\;1)", R"(\tau )", R"(\frac{1}{2},)", R"(\sum_{i = 0}^{\infty}\pi)"}) EXPECT_FALSE(isGlyphString(c));
}

TEST_F(ParseTeXTest, isTeXCommandEnd)
{
   std::string cmd;

   // Test single non-letter character commands
   for(const std::string cmd0: {"\\;", "\\\\", "\\{", "\\[", "\\%"})
      for(const auto c: {' ', ',', ';', ':', '(', '\t', '\n', '\\', '1', 'a'})
      {
         cmd = cmd0 + c + ' ';
         EXPECT_TRUE(isTeXCommandEnd(cmd.end() - 2, cmd.end()));
      }

   // Test single letter character and word commands
   for(const std::string cmd0: {"\\p", "\\tau"})
   {
      for(const auto c: {' ', ',', ';', ':', '(', '\t', '\n', '\\', '1'})
      {
         cmd = cmd0 + c + ' ';
         EXPECT_TRUE(isTeXCommandEnd(cmd.end() - 2, cmd.end()));
      }
      cmd = cmd0 + 'a' + ' ';
      EXPECT_FALSE(isTeXCommandEnd(cmd.end() - 2, cmd.end()));
   }
}

TEST_F(ParseTeXTest, isTeXCharCommand)
{
   for(const std::string cmd: {"\\ ", "\\p", "\\!", "\\, ", "\\%,", "\\$:", "\\;end", "\\ end"})
   {
      const auto [is_char, end_iter] = isTeXCharCommand(cmd.begin(), cmd.end());
      EXPECT_TRUE(is_char);
      EXPECT_EQ(end_iter, cmd.begin() + 2);
   }
   for(const std::string cmd: {"\\", "\\pi", "\\p_1", "\\frac{1}{2}"})
   {
      const auto [is_char, end_iter] = isTeXCharCommand(cmd.begin(), cmd.end());
      EXPECT_FALSE(is_char);
      EXPECT_EQ(end_iter, cmd.end());
   }
   for(const std::string cmd: {" ", "2", "pi", "p_1", "frac{1}{2}"}) EXPECT_DEATH(isTeXCharCommand(cmd.begin(), cmd.end()), "");
}

TEST_F(ParseTeXTest, isTeXWordCommand)
{
   for(const std::string cmd: {"\\pi", "\\tau_0", "\\frac{1}{2}", "\\bm"})
   {
      const auto [is_word, end_iter] = isTeXWordCommand(cmd.begin(), cmd.end());
      const auto end_iter_check = std::find_if_not(cmd.begin() + 1, cmd.end(), [](auto c) { return std::isalpha(c); });
      EXPECT_TRUE(is_word);
      EXPECT_EQ(end_iter, end_iter_check);
   }
   for(const std::string cmd: {"\\ ", "\\p", "\\!", "\\, ", "\\%,", "\\$:", "\\;end", "\\ end"})
   {
      const auto [is_word, end_iter] = isTeXWordCommand(cmd.begin(), cmd.end());
      EXPECT_FALSE(is_word);
      EXPECT_EQ(end_iter, cmd.end());
   }
   for(const std::string cmd: {" ", "2", "pi", "p_1", "frac{1}{2}"}) EXPECT_DEATH(isTeXWordCommand(cmd.begin(), cmd.end()), "");
}

TEST_F(ParseTeXTest, GetTeXCommandInfo)
{
   // Single-character TeX commands
   for(const std::string cmd: {"\\%", "\\& Sons", "\\#123", "\\^2", "\\_i", "\\@gmail.com"})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      EXPECT_EQ(cmd_prefix_end, cmd.begin() + 2);
      EXPECT_FALSE(trailing_args);
      EXPECT_TRUE(not_spacing);
   }
   for(const std::string cmd: {"\\!", "\\, ", "\\: 5", "\\> x^2", "\\; (x + 1)", "\\ end."})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      EXPECT_EQ(cmd_prefix_end, cmd.begin() + 2);
      EXPECT_FALSE(trailing_args);
      EXPECT_FALSE(not_spacing);
   }

   // Word TeX commands with no trailing arguments
   auto find_end_iter = [](const std::string& cmd){ return std::find_if_not(cmd.begin() + 1, cmd.end(), [](auto c){return std::isalpha(c);}); };
   for(const std::string cmd: {"\\pi", "\\tau ", "\\gamma, ", "\\beta.", "\\eta: 5x", "\\epsilon = 1e-5"})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      const auto end_iter = find_end_iter(cmd);

      EXPECT_EQ(cmd_prefix_end, end_iter);
      EXPECT_FALSE(trailing_args);
      EXPECT_TRUE(not_spacing);
   }
   for(const std::string cmd: {"\\quad", "\\qquad", "\\quad,", "\\qquad ", "\\quad5x", "\\qquad x = 1."})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      const auto end_iter = find_end_iter(cmd);

      EXPECT_EQ(cmd_prefix_end, end_iter);
      EXPECT_FALSE(trailing_args);
      EXPECT_FALSE(not_spacing);
   }

   // Word TeX commands with trailing arguments/script text
   for(const std::string cmd: {"\\frac{1}{2}", "\\dotp{x}{y}", "\\tau_0 ", "\\sum_{i=0} x", "\\intt_{\\Gamma} x dx"})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      const auto end_iter = find_end_iter(cmd);

      EXPECT_EQ(cmd_prefix_end, end_iter);
      EXPECT_TRUE(trailing_args);
      EXPECT_TRUE(not_spacing);
   }
   for(const std::string cmd: {"\\dec{1}" , "\\decc{1} " , "\\deccc{1}, " , "\\inc{1}: " , "\\incc{1};" , "\\inccc{1}(x)",
                               "\\vdec{1}", "\\vdecc{1}.", "\\vdeccc{1} ", "\\vinc{1} = 5", "\\vincc{1} ", "\\vinccc{1}"})
   {
      const auto [cmd_prefix_end, trailing_args, not_spacing] = GetTeXCommandInfo(cmd.begin(), cmd.end());
      const auto end_iter = find_end_iter(cmd);

      EXPECT_EQ(cmd_prefix_end, end_iter);
      EXPECT_TRUE(trailing_args);
      EXPECT_FALSE(not_spacing);
   }

   // Death tests
   for(const std::string cmd: {"x", "1", "\\1", "\\"}) EXPECT_DEATH(GetTeXCommandInfo(cmd.begin(), cmd.end()), "");
}

TEST_F(ParseTeXTest, ParseTeXScriptText)
{
   for(const std::string str: {"x^2", "\\tau_0 ", "(x + 1)^3 + 1", "\\sum_i x", "\\sum_{i} n", "c_{5}x^5", "e^{x}"})
   {
      const auto iter_start = std::find_if(str.begin(), str.end(), [](auto c){return c == OneOf('_', '^');});
      auto iter = iter_start;
      const auto objects = ParseTeXScriptText(iter, str.end());
      EXPECT_EQ(objects.size(), 1);

      const auto glyph = std::dynamic_pointer_cast<TeXGlyph>(objects[0]);
      EXPECT_TRUE(glyph->Rendered());

      if(*(iter_start + 1) != '{') // Not braced
      {
         EXPECT_EQ(glyph->Text(), std::string(1, *(iter_start + 1)));
         EXPECT_EQ(iter, iter_start + 2);
      }
      else // Braced
      {
         EXPECT_EQ(glyph->Text(), std::string(1, *(iter_start + 2)));
         EXPECT_EQ(iter, iter_start + 4);
      }
   }
}

}

#endif