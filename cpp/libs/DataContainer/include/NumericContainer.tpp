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

namespace aprn {
namespace detail {

/** Scalar arithmetic operator overloads. */
template<Arithmetic T, class D>
constexpr D
NumericContainer<T, D>::operator+(const std::convertible_to<T> auto scalar) const
{
  D out(Derived());
  FOR_EACH(entry, out) entry += scalar;
  return out;
}

template<Arithmetic T, class D>
constexpr D
NumericContainer<T, D>::operator-(const std::convertible_to<T> auto scalar) const
{
  D out(Derived());
  FOR_EACH(entry, out) entry -= scalar;
  return out;
}

template<Arithmetic T, class D>
constexpr D
NumericContainer<T, D>::operator*(const std::convertible_to<T> auto scalar) const
{
  D out(Derived());
  FOR_EACH(entry, out) entry *= scalar;
  return out;
}

template<Arithmetic T, class D>
constexpr D
NumericContainer<T, D>::operator/(const std::convertible_to<T> auto scalar) const
{
  DEBUG_ASSERT(!isEqual(scalar, Zero), "Cannot divide by zero.")
  D out(Derived());
  FOR_EACH(entry, out) entry /= scalar;
  return out;
}

template<Arithmetic T, class D>
constexpr D&
NumericContainer<T, D>::operator+=(const std::convertible_to<T> auto scalar)
{
  FOR_EACH(entry, Derived()) entry += scalar;
  return Derived();
}

template<Arithmetic T, class D>
constexpr D&
NumericContainer<T, D>::operator-=(const std::convertible_to<T> auto scalar)
{
  FOR_EACH(entry, Derived()) entry -= scalar;
  return Derived();
}

template<Arithmetic T, class D>
constexpr D&
NumericContainer<T, D>::operator*=(const std::convertible_to<T> auto scalar)
{
  FOR_EACH(entry, Derived()) entry *= scalar;
  return Derived();
}

template<Arithmetic T, class D>
constexpr D&
NumericContainer<T, D>::operator/=(const std::convertible_to<T> auto scalar)
{
  DEBUG_ASSERT(!isEqual(scalar, Zero), "Cannot divide by zero.")
  FOR_EACH(entry, Derived()) entry /= scalar;
  return Derived();
}

/** Entry-wise binary arithmetic operator overloads. */
template<Arithmetic T, class D>
template<class D2>
constexpr D
NumericContainer<T, D>::operator+(const NumericContainer<T, D2>& container) const
{
  D out(Derived());
  FOR(i, out.size()) out[i] += container.Derived()[i];
  return out;
}

template<Arithmetic T, class D>
template<class D2>
constexpr D
NumericContainer<T, D>::operator-(const NumericContainer<T, D2>& container) const
{
  D out(Derived());
  FOR(i, out.size()) out[i] -= container.Derived()[i];
  return out;
}

template<Arithmetic T, class D>
template<class D2>
constexpr D
NumericContainer<T, D>::operator*(const NumericContainer<T, D2>& container) const
{
  D out(Derived());
  FOR(i, out.size()) out[i] *= container.Derived()[i];
  return out;
}

template<Arithmetic T, class D>
template<class D2>
constexpr D
NumericContainer<T, D>::operator/(const NumericContainer<T, D2>& container) const
{
  D out(Derived());

  FOR(i, out.size())
  {
    DEBUG_ASSERT(!isEqual(container.Derived()[i], Zero), "Cannot divide by zero.")
    out[i] /= container.Derived()[i];
  }
  return out;
}

template<Arithmetic T, class D>
template<class D2>
constexpr D&
NumericContainer<T, D>::operator+=(const NumericContainer<T, D2>& container)
{
  FOR(i, Derived().size()) Derived()[i] += container.Derived()[i];
  return Derived();
}

template<Arithmetic T, class D>
template<class D2>
constexpr D&
NumericContainer<T, D>::operator-=(const NumericContainer<T, D2>& container)
{
  FOR(i, Derived().size()) Derived()[i] -= container.Derived()[i];
  return Derived();
}

template<Arithmetic T, class D>
template<class D2>
constexpr D&
NumericContainer<T, D>::operator*=(const NumericContainer<T, D2>& container)
{
  FOR(i, Derived().size()) Derived()[i] *= container.Derived()[i];
  return Derived();
}

template<Arithmetic T, class D>
template<class D2>
constexpr D&
NumericContainer<T, D>::operator/=(const NumericContainer<T, D2>& container)
{
  FOR(i, Derived().size())
  {
    DEBUG_ASSERT(!isEqual(container.Derived()[i], Zero), "Cannot divide by zero.")
    Derived()[i] /= container.Derived()[i];
  }
  return Derived();
}

/** Entry-wise unary operator overloads. */
template<Arithmetic T, class D>
constexpr D
NumericContainer<T, D>::operator-() const { return -1 * Derived(); }

/** Entry randomisation. */
template<Arithmetic T, class D>
void NumericContainer<T, D>::Randomise() { FOR_EACH(entry, Derived()) entry = Randomiser(); }

template<Arithmetic T, class D>
void NumericContainer<T, D>::ResetRandomiser(const T min, const T max) { Randomiser.Reset(min, max); }

}//Detail

/** Stand-alone Operator overloads. */
template<Arithmetic T, class D>
constexpr D
operator*(const std::convertible_to<T> auto scalar, const detail::NumericContainer<T, D>& container) { return container.Derived() * scalar; }

}//aprn