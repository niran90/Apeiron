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

#include "Constants.h"
#include "Debug.h"
#include "Types.h"

#include <chrono>
#include <random>

namespace aprn{

/***************************************************************************************************************************************************************
* Abstract Random Number Base Class
***************************************************************************************************************************************************************/
template<typename T, TypeCategory cat = GetTypeCategory<T>()>
class RandomBase
{
protected:
  RandomBase() : Generator_(Device_()) {}

public:
  /** Subscript operator overload to generate a random number. */
  virtual T operator()() = 0;

protected:
  std::random_device Device_;
  std::mt19937       Generator_;
};

/** Forward declaration of Random class and template deduction guide. */
template<typename T, TypeCategory cat = GetTypeCategory<T>()> class Random;
template<typename T> Random(T, T) -> Random<T>;

/***************************************************************************************************************************************************************
* Random Boolean Class
***************************************************************************************************************************************************************/
template<typename T>
class Random<T, TypeCategory::Boolean> : public RandomBase<T, TypeCategory::Boolean>
{
public:
  Random() : Distribution_(0, 1) {}

  /** Overloaded subscript operator to generate a random boolean. */
  inline T operator()() override { return Distribution_(this->Generator_); }

private:
  std::uniform_int_distribution<UChar> Distribution_;
};

/***************************************************************************************************************************************************************
* Random Integer Class
***************************************************************************************************************************************************************/
template<typename T>
class Random<T, TypeCategory::Integral> : public RandomBase<T, TypeCategory::Integral>
{
public:
  Random() : Distribution_(-1, 1) {}

  Random(const T min, const T max) : Distribution_(min, max) {}

  /** Overloaded subscript operator to generate a random integer. */
  inline T operator()() override { return Distribution_(this->Generator_); }

  /** Reset min/max bounds for the distribution. */
  void Reset(const T min, const T max) { Distribution_ = std::uniform_int_distribution<T>(min, max); }

private:
  std::uniform_int_distribution<T> Distribution_;
};

/***************************************************************************************************************************************************************
* Random Floating-point Class
***************************************************************************************************************************************************************/
template<typename T>
class Random<T, TypeCategory::FloatingPoint> : public RandomBase<T, TypeCategory::FloatingPoint>
{
public:
  Random() : Distribution_(-One, One) {}

  Random(const T min, const T max) : Distribution_(min, max) {}

  /** Overloaded subscript operator to generate a random floating-point value. */
  inline T operator()() override { return Distribution_(this->Generator_); }

  /** Reset min/max bounds for the distribution. */
  void Reset(const T min, const T max) { Distribution_ = std::uniform_real_distribution<T>(min, max); }

private:
  std::uniform_real_distribution<T> Distribution_;
};

}
