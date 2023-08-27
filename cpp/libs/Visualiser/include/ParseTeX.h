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
#include "TeXGlyph.h"
#include "TeXObject.h"

#include <string>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* TeX Parsing Functions
***************************************************************************************************************************************************************/
SPtr<TeXGlyph> ParseTeXChar(const char c);

DArray<SPtr<TeXObject>> ParseTeXString(const std::string& str);

template<CharIterator It>
DArray<SPtr<TeXObject>> ParseTeXString(It& current, const It last, const bool is_math_mode = false);

template<CharIterator It>
SPtr<TeXObject> ParseTeXObject(It& current, const It last, const bool is_math_mode = false);

template<CharIterator It>
SPtr<TeXObject> ParseTeXCommand(It& current, const It last, const bool is_math_mode = false);

template<CharIterator It>
DArray<SPtr<TeXObject>> ParseTeXMath(It& current, const It last);

template<CharIterator It>
std::pair<DArray<SPtr<TeXObject>>, std::string> ParseAllTeXScriptText(It& current, const It last);

/***************************************************************************************************************************************************************
* TeX Parsing Helper Functions
***************************************************************************************************************************************************************/
bool isGlyphString(const std::string_view& tex_str);

template<CharIterator It>
bool isTeXCommandEnd(const It current, const It last);

template<CharIterator It>
std::pair<bool, It> isTeXCharCommand(const It first, const It last);

template<CharIterator It>
std::pair<bool, It> isTeXWordCommand(const It first, const It last);

template<CharIterator It>
std::tuple<It, bool, bool> GetTeXCommandInfo(const It current, const It last);

template<CharIterator It>
DArray<SPtr<TeXObject>> ParseTeXScriptText(It& current, const It last);

template<CharIterator It>
It GetTeXCommandPrefixEnd(const It first, const It last, const bool is_char_cmd);

}

#include "ParseTeX.tpp"
