#pragma once

#include "../../../include/Global.h"

#include <array>
#include <initializer_list>
#include <vector>

namespace Apeiron {

/***************************************************************************************************************************************************************
* Array Abstract Base Class
***************************************************************************************************************************************************************/
template <typename T, class derived>
class Array
{
protected:
  constexpr Array() = default;

public:
  /** Size and Index Range-checking */
  constexpr void IndexBoundCheck(const std::size_t _index) const;

  constexpr void SizeCheck(const std::size_t _size0, const std::size_t _size1) const;

  /** Subscript Operator Overloads */
  constexpr T& operator[](const std::size_t _index);

  constexpr const T& operator[](const std::size_t _index) const;

  /** Assignment Operator Overloads */
  template <class t_rhs_type>
  constexpr derived& operator=(const t_rhs_type& _value) noexcept;

  constexpr derived& operator=(const std::initializer_list<T>& _value_list) noexcept;

private:
  /** Derived Class Access */
  constexpr derived& Derived() noexcept { return static_cast<derived&>(*this); }

  constexpr const derived& Derived() const noexcept { return static_cast<const derived&>(*this); }
};

/** Non-member functions */
template <typename T, class derived>
std::ostream& operator<<(std::ostream& _output_stream, const Array<T, derived>& _array_base);

/***************************************************************************************************************************************************************
* Static Array Class
***************************************************************************************************************************************************************/
template <typename T, std::size_t N>
class StaticArray : public std::array<T, N>,
                    public Array<T, StaticArray<T, N>>
{
  using Base = Array<T, StaticArray<T, N>>;

public:
  /** Constructors */
  constexpr StaticArray();

  constexpr StaticArray(const T& _value);

  constexpr StaticArray(const std::initializer_list<T>& _list);

  template <class iter>
  constexpr StaticArray(const iter _first, const iter _last);

  /** Operators */
  using Base::operator[];
  using Base::operator=;

  friend Base;
};

/***************************************************************************************************************************************************************
* Dynamic Array Class
***************************************************************************************************************************************************************/
template <typename T>
class DynamicArray : public std::vector<T>,
                     public Array<T, DynamicArray<T>>
{
  using Base = Array<T, DynamicArray<T>>;

public:
  /** Constructors */
  DynamicArray();

  DynamicArray(const std::size_t _size);

  DynamicArray(const std::size_t _size, const T& _value);

  DynamicArray(const std::initializer_list<T>& _list);

  template <class iter>
  DynamicArray(const iter _first, const iter _last);

  /** Operators */
  using Base::operator[];
  using Base::operator=;

  friend Base;
};

}

#include "../src/Array.cpp"
