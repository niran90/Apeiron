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

namespace aprn::flmgr {

template<bool wide>
BaseFile<wide>::BaseFile(const Path& file_path, ModeType auto... open_modes)
   : Path_(file_path), Stream_(file_path,
                               (static_cast<std::ios_base::openmode>(open_modes) | ...)), Modes_({open_modes...}) { Init(file_path); }

template<bool wide>
void
BaseFile<wide>::Open(const Path& file_path, ModeType auto... open_modes)
{
   ASSERT(!isOpen(), "The file ", file_path.filename(), " is already open.")

   Path_  = file_path;
   Modes_ = {open_modes... };
   Stream_.open(file_path, (static_cast<std::ios_base::openmode>(open_modes) | ...));

   Init(file_path);
}

template<bool wide>
template<char sep, typename T, typename... Ts>
void
BaseFile<wide>::Read(T& data, Ts&... trailing_data)
{
   DEBUG_ASSERT(isReadable(), "The file must be readable to read data.")

   Stream_ >> data;
   if constexpr(sizeof...(Ts))
   {
      if constexpr(sep != '\0')
      {
         DEBUG_ASSERT(Stream_.peek() == sep, "The data to be read in is not delimitted by a '", sep, "' character.")
         Stream_ >> Unmove(char());
      }
      Read<sep>(trailing_data...);
   }
}

template<bool wide>
template<char sep, typename T, typename... Ts>
void
BaseFile<wide>::Write(const T& data, const Ts&... trailing_data)
{
   DEBUG_ASSERT(isWritable(), "The file must be writable to write a line.")

   Stream_ << data;
   if constexpr(sizeof...(Ts))
   {
      if constexpr(sep != '\0') Stream_ << sep;
      Write<sep>(trailing_data...);
   }
}

}