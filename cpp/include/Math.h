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

#include "Comparators.h"
#include "Constants.h"
#include "Debug.h"
#include "Types.h"

#include <cmath>
#include <functional>
#include <numeric>

namespace aprn{

/***************************************************************************************************************************************************************
* Math Support Functions
***************************************************************************************************************************************************************/
namespace detail {

/** Modulo function for integer types. */
template<typename T>
constexpr T
Modulo(const T numerator, const T denominator, std::true_type, std::false_type) { return numerator % denominator; }

/** Modulo function for floating-point types. */
template<typename T>
constexpr T
Modulo(const T numerator, const T denominator, std::false_type, std::true_type) { return std::fmod(numerator, denominator); }

}

/***************************************************************************************************************************************************************
* Sequences and Series
***************************************************************************************************************************************************************/
template<typename T, T... sequence>
constexpr T
GetInteger(std::integer_sequence<T, sequence...>, size_t index)
{
  constexpr T sequence_array[] = {sequence...};
  return sequence_array[index];
}

/***************************************************************************************************************************************************************
* General Arithmetic Functions
***************************************************************************************************************************************************************/
/** Sum the terms of a sequence together. */
template<typename... T>
constexpr auto
Sum(const T... values) { return (values + ... + 0); }

/** Sum the terms of a sequence between two iterators. */
template<class It>
constexpr auto
Sum(const It first, const It last)
{
  typedef typename std::iterator_traits<It>::value_type data_type;
  return std::accumulate(first, last, static_cast<data_type>(0));
}

/** Product the terms of a sequence with each other. */
template<typename... T>
constexpr auto
Product(const T... values) { return (values * ... * 1); }

/** Product the terms of a sequence between two iterators. */
template<class It>
constexpr auto
Product(const It first, const It last)
{
  typedef typename std::iterator_traits<It>::value_type data_type;
  return std::accumulate(first, last, static_cast<data_type>(1), std::multiplies<data_type>());
}

/** Division function. */
template<typename T>
constexpr Real
Divide(const T _numerator, const T _denominator)
{
  return !isEqual(static_cast<Real>(_denominator), Zero) ? static_cast<Real>(_numerator) / static_cast<Real>(_denominator) :
         throw std::logic_error("Denominator must be non-zero during division.");
}

/** Modulo function. */
template<typename T>
constexpr T
Modulo(const T _numerator, const T _denominator)
{
  return !isEqual(static_cast<Real>(_denominator), Zero) ?
         detail::Modulo<T>(_numerator, _denominator, std::is_integral<T>(), std::is_floating_point<T>()) :
         throw std::logic_error("Denominator must be non-zero during modulo operation.");
}

/***************************************************************************************************************************************************************
* Combinatorial Functions
***************************************************************************************************************************************************************/
/** Factorial function (i.e. n!) */
constexpr unsigned int
Factorial(const unsigned int _x)
{
  return _x <= 20 ? (0 < _x ? _x*Factorial(_x - 1) : 1) :
         throw std::logic_error("Cannot currently compute the factorial of a number larger than 20.");
}

/** Factorial quotient function (i.e. n!/m!) */
constexpr unsigned int
FactorialQuotient(const unsigned int _numerator, const unsigned int _denominator)
{
  return _numerator <= 20 && _denominator <= 20 ?
         (_denominator <= _numerator ? (_denominator < _numerator ? _numerator * FactorialQuotient(_numerator - 1, _denominator) : 1) :
         throw std::logic_error("Numerator must be larger than the denominator.")) :
         throw std::logic_error("Cannot currently compute the factorial of a number larger than 20.");
}

/** Compination function ("n choose r"). */
constexpr unsigned int
Choose(const unsigned int _n, const unsigned int _r) {
  return _r <= _n ? FactorialQuotient(_n, _n - _r) / Factorial(_r) : throw std::logic_error("Numerator must be larger than the denominator.");
}

/***************************************************************************************************************************************************************
* Power Functions
***************************************************************************************************************************************************************/
/** Constexpr power function with integral exponents [significantly faster than std::pow(x, i)]. */
template<typename T>
constexpr T
iPow(const T x, const unsigned int exponent)
{
  return exponent <= 30 ? (exponent == 0 ? static_cast<T>(1) : x * iPow(x, exponent - 1)) :
         throw std::logic_error("Cannot currently compute the power with an exponent larger than 30.");
}

/** Square of a value [significantly faster than std::pow(x, 2)]. */
template<typename T>
constexpr T
Square(const T x) { return iPow(x, 2); }

/** Cube of a value [significantly faster than std::pow(x, 3)]. */
template<typename T>
constexpr T
Cube(const T x) { return iPow(x, 3); }

}//aprn

