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
#include "../../../include/Concepts.h"
#include "../../../include/Random.h"

namespace aprn {
namespace detail {

/***************************************************************************************************************************************************************
* Numeric Data Container Class and Arithmetic Operations
***************************************************************************************************************************************************************/
template<Arithmetic T, class D>
class NumericContainer
{
protected:
  constexpr NumericContainer() {}

public:
  /** Scalar arithmetic operator overloads. */
  constexpr D operator+(const std::convertible_to<T> auto scalar) const;

  constexpr D operator-(const std::convertible_to<T> auto scalar) const;

  constexpr D operator*(const std::convertible_to<T> auto scalar) const;

  constexpr D operator/(const std::convertible_to<T> auto scalar) const;

  constexpr D& operator+=(const std::convertible_to<T> auto scalar);

  constexpr D& operator-=(const std::convertible_to<T> auto scalar);

  constexpr D& operator*=(const std::convertible_to<T> auto scalar);

  constexpr D& operator/=(const std::convertible_to<T> auto scalar);

  /** Entry-wise binary arithmetic operator overloads. */
  template<class D2>
  constexpr D operator+(const NumericContainer<T, D2>& container) const;

  template<class D2>
  constexpr D operator-(const NumericContainer<T, D2>& container) const;

  template<class D2>
  constexpr D operator*(const NumericContainer<T, D2>& container) const;

  template<class D2>
  constexpr D operator/(const NumericContainer<T, D2>& container) const;

  template<class D2>
  constexpr D& operator+=(const NumericContainer<T, D2>& container);

  template<class D2>
  constexpr D& operator-=(const NumericContainer<T, D2>& container);

  template<class D2>
  constexpr D& operator*=(const NumericContainer<T, D2>& container);

  template<class D2>
  constexpr D& operator/=(const NumericContainer<T, D2>& container);

  /** Entry-wise unary operator overloads. */
  constexpr D operator-() const;

  /** Entry randomisation. */
  void Randomise();

  static void ResetRandomiser(const T min, const T max);

  /** Derived class access. */
  constexpr D&
  Derived() noexcept { return static_cast<D&>(*this); }

  constexpr const D&
  Derived() const noexcept { return static_cast<const D&>(*this); }

private:
  inline static Random<T> Randomiser = Random<T>();
};

}//detail

/** Stand-alone Operator overloads. */
template<Arithmetic T, class D>
constexpr D operator*(const std::convertible_to<T> auto scalar, const detail::NumericContainer<T, D>& container);

}//aprn

#include "NumericContainer.tpp"
