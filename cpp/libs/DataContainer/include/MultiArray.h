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
#include "Array.h"

namespace aprn{

/***************************************************************************************************************************************************************
* Multi-dimensional Array Abstract Base Class
***************************************************************************************************************************************************************/
template<typename T, class D>
class MultiArray
{
protected:
  constexpr MultiArray() = default;

public:
  /** Size and index range-checking. */
  constexpr void
  MultiIndexBoundCheck(const std::convertible_to<size_t> auto... multi_index) const;

  /** Multi-dimensional subscript index toggling. */
  constexpr size_t
  ComputeLinearIndex(const std::convertible_to<size_t> auto... multi_index) const;

  constexpr auto
  ComputeMultiIndex(size_t index) const;

  /** Operator overloads. */
  constexpr T&
  operator()(std::convertible_to<size_t> auto... multi_index);

  constexpr const T&
  operator()(const std::convertible_to<size_t> auto... multi_index) const;

  constexpr D&
  operator=(const std::initializer_list<T>& _value_array) noexcept;

  constexpr D&
  operator=(const std::initializer_list<std::initializer_list<T>>& _value_matrix) noexcept;

  /** Iterators */
  constexpr auto
  begin() { return Derived().Entries.begin(); }

  constexpr auto
  begin() const { return Derived().Entries.begin(); }

  constexpr auto
  end() { return Derived().Entries.end(); }

  constexpr auto
  end() const { return Derived().Entries.end(); }

private:
  /** Derived class access. */
  constexpr D&
  Derived() noexcept { return static_cast<D&>(*this); }

  constexpr const D&
  Derived() const noexcept { return static_cast<const D&>(*this); }
};

/***************************************************************************************************************************************************************
* Static Multi-dimensional Array Class
***************************************************************************************************************************************************************/
template<typename T, size_t ...dims>
class StaticMultiArray : public MultiArray<T, StaticMultiArray<T, dims...>>
{
public:
  /** Constructors. */
  constexpr StaticMultiArray();

  explicit constexpr StaticMultiArray(const T value);

private:
  constexpr static StaticArray<size_t, sizeof...(dims)> Dimensions{dims...};
  constexpr static size_t nEntries{Product(dims...)};
  StaticArray<T, nEntries> Entries;

  friend MultiArray<T, StaticMultiArray<T, dims...>>;
};

/***************************************************************************************************************************************************************
* Dynamic Multi-dimensional Array Class
***************************************************************************************************************************************************************/
template<typename T>
class DynamicMultiArray : public MultiArray<T, DynamicMultiArray<T>>
{
public:
  /** Constructors. */
  DynamicMultiArray();

  DynamicMultiArray(const std::convertible_to<size_t> auto... _dimensions);

  /** Multi-array resize. */
  void Resize(const std::convertible_to<size_t> auto... _dimensions);

private:
  DynamicArray<size_t> Dimensions;
  size_t nEntries;
  DynamicArray<T> Entries;

  friend MultiArray<T, DynamicMultiArray<T>>;
};

}

#include "MultiArray.tpp"
