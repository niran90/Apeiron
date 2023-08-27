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

namespace aprn {

/***************************************************************************************************************************************************************
* Multi-dimensional Array Abstract Base Class
***************************************************************************************************************************************************************/

/** Size and index range-checking */
template<typename T, class D>
constexpr void
MultiArray<T, D>::MultiIndexBoundCheck(const std::convertible_to<size_t> auto... multi_index) const
{
#ifdef DEBUG_MODE
  ASSERT(sizeof...(multi_index) == Derived().Dimensions.size(), "Multi-index size mismatch.")

  auto& dims = Derived().Dimensions;
  const size_t indices[] = {static_cast<size_t>(multi_index)...};
  FOR(i, dims.size()) ASSERT(indices[i] < dims[i], "Multi index component ", indices[i], " must be lesser than ", dims[i], ".")
#endif
}

/** Multi-dimensional subscript index toggling */
template<typename T, class D>
constexpr size_t
MultiArray<T, D>::ComputeLinearIndex(const std::convertible_to<size_t> auto... multi_index) const
{
  MultiIndexBoundCheck(multi_index...);

  auto& dims = Derived().Dimensions;
  const size_t indices[] = {static_cast<size_t>(multi_index)...};

  size_t index(0), factor(1);
  FOR(i, dims.size())
  {
    index += factor * indices[i];
    factor *= dims[i];
  }
//  this->IndexBoundCheck(index);

  return index;
}

template<typename T, class D>
constexpr auto
MultiArray<T, D>::ComputeMultiIndex(size_t index) const
{
  EXIT("Yet to complete.")

//  constexpr auto& dims = Derived().Dimensions;
//  constexpr auto& n_entries = Derived().nEntries;
//
//  DEBUG_ASSERT(isBounded(index, 0, n_entries), "The passed index ", index, " should be in the range [0, ", n_entries, "].")
//
//  StaticArray<size_t, dims.size()> res;

//  size_t factor = Multiply(dims...);
//  for (size_t i = sizeof...(dims); i != 0; --i) {
//    factor /= dims[i - 1];
//    res[i - 1] = index / factor;
//    assert(res[i - 1] < dims[i - 1]);
//    index -= res[i - 1] * factor;
//  }
//  return res;
  return 0;
}

/** Operator overloads. */
template<typename T, class D>
constexpr T&
MultiArray<T, D>::operator()(std::convertible_to<size_t> auto... multi_index)
{
  return Derived().Entries[ComputeLinearIndex(multi_index...)];
}

template<typename T, class D>
constexpr const T&
MultiArray<T, D>::operator()(const std::convertible_to<size_t> auto... multi_index) const
{
  return Derived().Entries[ComputeLinearIndex(multi_index...)];
}

template<typename T, class D>
constexpr D&
MultiArray<T, D>::operator=(const std::initializer_list<T>& _value_array) noexcept
{
  D& derived_class = Derived();

  SizeCheck(derived_class.Dimensions.size(), 1);
  SizeCheck(_value_array.size(), derived_class.Dimensions[0].size());
  FOR(i, _value_array.size()) derived_class.Entries[i] = _value_array[i];

  return derived_class;
}

template<typename T, class D>
constexpr D&
MultiArray<T, D>::operator=(const std::initializer_list<std::initializer_list<T>>& _value_matrix) noexcept
{
  D& derived_class = Derived();

  SizeCheck(derived_class.Dimensions.size(), 2);
  SizeCheck(_value_matrix.size(), derived_class.Dimensions[0].size());

  FOR(i, _value_matrix.size())
  {
    SizeCheck(_value_matrix[i].size(), derived_class.Dimensions[1].size());
    FOR(j, _value_matrix[i].size()) derived_class(i, j) = _value_matrix[i][j];
  }

  return derived_class;
}

/***************************************************************************************************************************************************************
* Static Multi-dimensional Array Class
***************************************************************************************************************************************************************/
template<typename T, size_t... dims>
constexpr StaticMultiArray<T, dims...>::StaticMultiArray()
  : StaticMultiArray(StaticInitValue<T>()) {}

template<typename T, size_t ...dims>
constexpr StaticMultiArray<T, dims...>::StaticMultiArray(const T value)
  : Entries(value)
{
  STATIC_ASSERT(0 < sizeof...(dims), "A multi-dimensional array must have at least 1 dimension.")
}

/***************************************************************************************************************************************************************
* Dynamic Multi-dimensional Array Class
***************************************************************************************************************************************************************/

/** Constructors/Destructors */
template<typename T>
DynamicMultiArray<T>::DynamicMultiArray()
  : DynamicMultiArray(0) {}

template<typename T>
DynamicMultiArray<T>::DynamicMultiArray(const std::convertible_to<size_t> auto... _dimensions)
  : Dimensions{static_cast<size_t>(_dimensions)...}, nEntries(Product(_dimensions...)), Entries(nEntries, DynamicInitValue<T>()) {}

/** Multi-array Resize Functions */
template<typename T>
void DynamicMultiArray<T>::Resize(const std::convertible_to<size_t> auto... _dimensions)
{
  Dimensions = {_dimensions...};
  nEntries = Product(_dimensions...);
  Entries.resize(nEntries, DynamicInitValue<T>());
}

}
