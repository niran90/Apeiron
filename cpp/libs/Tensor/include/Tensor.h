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
#include "../../DataContainer/include/MultiArray.h"
#include "../../DataContainer/include/NumericContainer.h"

namespace aprn{

/***************************************************************************************************************************************************************
* Tensor Abstract Base Class
***************************************************************************************************************************************************************/
template<typename T, class D>
class Tensor : public detail::NumericContainer<T, Tensor<T, D>>
{
protected:
  constexpr Tensor();

public:
  /** Subscript operator overloads. */
  constexpr T&
  operator()(std::convertible_to<size_t> auto... multi_index);

  constexpr const T&
  operator()(const std::convertible_to<size_t> auto... multi_index) const;

  /** Assignment operator overloads. */
  constexpr D&
  operator=(const std::initializer_list<T>& _value_array) noexcept;

  constexpr D&
  operator=(const std::initializer_list<std::initializer_list<T>>& _value_matrix) noexcept;

  /** Iterators */
  constexpr auto
  begin() { return Derived().Entries.begin(); }

  constexpr const auto
  begin() const { return Derived().Entries.begin(); }

  constexpr auto
  end() { return Derived().Entries.end(); }

  constexpr const auto
  end() const { return Derived().Entries.end(); }

private:
  std::pair<size_t, size_t> Type;

  /** Derived class access. */
  constexpr D&
  Derived() noexcept { return static_cast<D&>(*this); }

  constexpr const D&
  Derived() const noexcept { return static_cast<const D&>(*this); }
};

/***************************************************************************************************************************************************************
* Static Tensor Class
***************************************************************************************************************************************************************/
template<typename T, size_t... dims>
class StaticTensor : public Tensor<T, StaticTensor<T, dims...>>
{
  friend Tensor<T, StaticTensor<T, dims...>>;

public:
  StaticTensor();

private:
  StaticMultiArray<T, dims...> Entries;
};

/***************************************************************************************************************************************************************
* Dynamic Tensor Class
***************************************************************************************************************************************************************/
template<typename T>
class DynamicTensor : public Tensor<T, DynamicTensor<T>>
{
  friend Tensor<T, DynamicTensor<T>>;

public:
  DynamicTensor();

  DynamicTensor(const std::convertible_to<size_t> auto... _dimensions);

  inline void Resize(const std::convertible_to<size_t> auto... _dimensions) { Entries.Resize(_dimensions...); }

private:
  DynamicMultiArray<T> Entries;
};


//StaticTensor<Real, 3> test;
//
//template<typename T, unsigned dim>
//class StaticVector : public StaticTensor<T, 1, dim> {};
//
//template<typename T, unsigned t_dimension0, unsigned t_dimension1>
//class StaticMatrix : public StaticTensor<T, 2, t_dimension0, t_dimension1> {};
//
//typedef StaticTensor<Real, 1, 1> Vector_StFl;

}

#include "../src/Tensor.cpp"
