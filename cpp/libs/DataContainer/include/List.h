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

#include "../include/Global.h"

#include <list>

namespace aprn {

template<typename T>
class List : public std::list<T>
{
public:
  List() = default;
  ~List() = default;

  /** Size and Index Range-checking ***************************************************************************************************************************/
  void IndexBoundCheck(const size_t index) const
  {
    DEBUG_ASSERT(this->size(), "The list has not yet been sized.")
    DEBUG_ASSERT(isBounded(index, size_t(0), this->size()), "The list index ", index, " must be in the range [0, ", this->size() - 1, "].")
  }

  void SizeCheck(const size_t _size0, const size_t _size1) const
  {
    DEBUG_ASSERT(areSizesEqual(_size0, _size1), "The list sizes ", _size0, " and ", _size1, " must be equal.")
  }

  /** Subscript Operator Overloads ****************************************************************************************************************************/
  T& operator[](const size_t index)
  {
    IndexBoundCheck(index);
    auto iterator = this->begin();
    FOR(i, index) ++iterator;
    return *iterator;
  }

  const T& operator[](const size_t index) const
  {
    IndexBoundCheck(index);
    auto iterator = this->begin();
    FOR(i, index) ++iterator;
    return *iterator;
  }

  /** Assignment Operator Overloads ***************************************************************************************************************************/
  List& operator=(const T& value) noexcept
  {
    FOR_EACH(entry, *this) entry = value;
    return *this;
  }

  List& operator=(const std::initializer_list<T>& value_list)
  {
    SizeCheck(value_list.size(), this->size());
    size_t index(0);
    FOR_EACH(entry, value_list) (*this)[index++] = entry;
    return *this;
  }
};

}
