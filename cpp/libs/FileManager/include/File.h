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
#include "FileSystem.h"

#include <fstream>
#include <optional>
#include <string>
#include <string_view>

namespace aprn::flmgr {

enum class Mode
{
   Read     = std::ios::in,     // Read
   Write    = std::ios::out,    // Write
   Append   = std::ios::app,    // Write
   Truncate = std::ios::trunc,  // Write
   AtEnd    = std::ios::ate,    // Read/Write
   Binary   = std::ios::binary, // Read/Write
};

template<typename T>
concept ModeType = isTypeSame<T, Mode>();

template<bool wide>
class BaseFile
{
   template<typename T1, typename T2> using ConditionalType = std::conditional_t<!wide, T1, T2>;

 public:
   BaseFile() = default;

   BaseFile(const Path& file_path, ModeType auto... open_modes);

   ~BaseFile();

   void Open(const Path& file_path, ModeType auto... open_modes);

   void Close();

   void Read(ConditionalType<std::string, std::wstring>& line);

   void Write(const ConditionalType<std::string_view, std::wstring_view>& str);

   template<char sep = '\0', typename T, typename... Ts>
   void Read(T& data, Ts&... trailing_data);

   template<char sep = '\0', typename T, typename... Ts>
   void Write(const T& data, const Ts&... trailing_data);

   inline bool isOpen() const { return Stream_.is_open(); }

   inline bool isValid() const { return Stream_.good(); }

   inline bool isEnd() const { return Stream_.eof(); }

   inline void SetPrecision(const unsigned n) { Stream_ << std::setprecision(n); }

 private:
   void Init(const Path& file_path);

   bool isReadable() const;

   bool isWritable() const;

   using FileStream = std::basic_fstream<ConditionalType<char, wchar_t>>;
   Path                 Path_;
   FileStream           Stream_;
   DArray<Mode>         Modes_;
   mutable Option<bool> Readable_;
   mutable Option<bool> Writable_;
};

typedef BaseFile<false> File;
typedef BaseFile<true>  WFile;

}

#include "File.tpp"
