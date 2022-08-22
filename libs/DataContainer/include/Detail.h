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

namespace aprn::detail {

/***************************************************************************************************************************************************************
* Static Array Support Functions
***************************************************************************************************************************************************************/

/** Initialise an std::array with a single value. */
template<typename T, size_t size>
constexpr auto InitStaticArray(const T& value)
{
  std::array<T, size> initialised_array;
  std::fill(initialised_array.begin(), initialised_array.end(), value);
  return initialised_array;
}

/** Initialise an std::array with an initializer list. */
template<typename T, size_t size>
constexpr auto InitStaticArray(const std::initializer_list<T>& list)
{
  std::array<T, size> initialised_array;
  std::copy(list.begin(), list.end(), initialised_array.begin());
  return initialised_array;
}

/** Initialise an std::array with a first and last iterator. */
template<typename T, size_t size, class iter>
constexpr auto InitStaticArray(const iter first, const iter last)
{
  std::array<T, size> initialised_array;
  std::copy(first, last, initialised_array.begin());
  return initialised_array;
}

}//aprn::detail