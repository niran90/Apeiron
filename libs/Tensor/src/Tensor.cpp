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

namespace aprn{

/***************************************************************************************************************************************************************
* Tensor Abstract Base Class
***************************************************************************************************************************************************************/
template<typename T, class derived>
constexpr Tensor<T, derived>::Tensor()
{

}

template<typename T, class derived>
constexpr T&
Tensor<T, derived>::operator()(std::convertible_to<size_t> auto... _multi_index)
{
  return Derived().Entries(_multi_index...);
}

template<typename T, class derived>
constexpr const T&
Tensor<T, derived>::operator()(const std::convertible_to<size_t> auto... _multi_index) const
{
  return Derived().Entries(_multi_index...);
}

template<typename T, class derived>
constexpr derived&
Tensor<T, derived>::operator=(const std::initializer_list<T>& _value_array) noexcept
{
  Derived().Entries = _value_array;
}

template<typename T, class derived>
constexpr derived&
Tensor<T, derived>::operator=(const std::initializer_list<std::initializer_list<T>>& _value_matrix) noexcept
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