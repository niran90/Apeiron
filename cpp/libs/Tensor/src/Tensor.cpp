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

#ifndef TENSOR_TEMPLATE_DEF
#define TENSOR_TEMPLATE_DEF

#include "../include/Tensor.h"

namespace aprn{

/***************************************************************************************************************************************************************
* Tensor Abstract Base Class
***************************************************************************************************************************************************************/
template<typename T, class D>
constexpr Tensor<T, D>::Tensor()
{

}

template<typename T, class D>
constexpr T&
Tensor<T, D>::operator()(std::convertible_to<size_t> auto... multi_index)
{
  return Derived().Entries(multi_index...);
}

template<typename T, class D>
constexpr const T&
Tensor<T, D>::operator()(const std::convertible_to<size_t> auto... multi_index) const
{
  return Derived().Entries(multi_index...);
}

template<typename T, class D>
constexpr D&
Tensor<T, D>::operator=(const std::initializer_list<T>& _value_array) noexcept
{
  Derived().Entries = _value_array;
}

template<typename T, class D>
constexpr D&
Tensor<T, D>::operator=(const std::initializer_list<std::initializer_list<T>>& _value_matrix) noexcept
{
  Derived().Entries = _value_matrix;
}

/***************************************************************************************************************************************************************
* Static Tensor Class
***************************************************************************************************************************************************************/
template<typename T, size_t... dims>
StaticTensor<T, dims...>::StaticTensor()
{
  STATIC_ASSERT(0 < sizeof...(dims), "A tensor must have at least 1 dimension.")
}

/***************************************************************************************************************************************************************
* Dynamic Tensor Class
***************************************************************************************************************************************************************/
template<typename T>
DynamicTensor<T>::DynamicTensor()
  : Entries() {}

template<typename T>
DynamicTensor<T>::DynamicTensor(const std::convertible_to<size_t> auto... _dimensions)
  : Entries(_dimensions...) {}

}

#endif