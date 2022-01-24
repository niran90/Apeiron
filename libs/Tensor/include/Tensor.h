#pragma once

#include "../../../include/Global.h"
#include "../../DataContainer/include/MultiArray.h"
#include "../../DataContainer/include/Detail.h"

namespace Apeiron{

/***************************************************************************************************************************************************************
* Tensor Abstract Base Class
***************************************************************************************************************************************************************/
template<typename T, class derived>
class Tensor : Detail::NumericContainer<T, Tensor<T, derived>>
{
protected:
  constexpr Tensor();

public:
  /** Subscript operator overloads. */
  constexpr T&
  operator()(std::convertible_to<std::size_t> auto ..._multi_index);

  constexpr const T&
  operator()(const std::convertible_to<std::size_t> auto ..._multi_index) const;

  /** Assignment operator overloads. */
  constexpr derived&
  operator=(const std::initializer_list<T>& _value_array) noexcept;

  constexpr derived&
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
  std::pair<std::size_t, std::size_t> Type;

  /** Derived class access. */
  constexpr derived&
  Derived() noexcept { return static_cast<derived&>(*this); }

  constexpr const derived&
  Derived() const noexcept { return static_cast<const derived&>(*this); }
};

/***************************************************************************************************************************************************************
* Static Tensor Class
***************************************************************************************************************************************************************/
template<typename T, std::size_t ...dimensions>
class StaticTensor : public Tensor<T, StaticTensor<T, dimensions...>>
{
  friend Tensor<T, StaticTensor<T, dimensions...>>;

public:
  StaticTensor();

private:
  StaticMultiArray<T, dimensions...> Entries;
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

  DynamicTensor(const std::convertible_to<std::size_t> auto ..._dimensions);

  inline void Resize(const std::convertible_to<std::size_t> auto ..._dimensions) { Entries.Resize(_dimensions...); }

private:
  DynamicMultiArray<T> Entries;
};


//StaticTensor<Float, 3> test;
//
//template<typename T, unsigned dim>
//class StaticVector : public StaticTensor<T, 1, dim> {};
//
//template<typename T, unsigned t_dimension0, unsigned t_dimension1>
//class StaticMatrix : public StaticTensor<T, 2, t_dimension0, t_dimension1> {};
//
//typedef StaticTensor<Float, 1, 1> Vector_StFl;

}

#include "../src/Tensor.cpp"
