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

#include <cmath>
#include <concepts>
#include <limits>
#include <memory>
#include <optional>
#include <type_traits>

namespace aprn {

/***************************************************************************************************************************************************************
* Type Categories
***************************************************************************************************************************************************************/
enum class TypeCategory
{
   Boolean,
   Integral,
   FloatingPoint,
   String,
   Other
};

/***************************************************************************************************************************************************************
* Boolean Type Wrapper (for std:: containers)
***************************************************************************************************************************************************************/
struct Bool
{
   constexpr Bool() : Bool(false) {}

   constexpr Bool(const bool value) : Value(value) {}

   constexpr operator bool() const { return Value; }

   friend std::ostream& operator<<(std::ostream& stream, const Bool value) { return stream << std::boolalpha << value.Value; }

 private:
   bool Value;
};

constexpr Bool True  = true;
constexpr Bool False = false;

/***************************************************************************************************************************************************************
* Integral Types
***************************************************************************************************************************************************************/

/** Char type definitions. */
typedef char                   Char;
typedef unsigned char          UChar;

/** Signed integer type definitions. */
typedef int                    Int;
typedef long int               LInt;
typedef long long int          LLInt;
typedef int8_t                 Int8;
typedef int16_t                Int16;
typedef int32_t                Int32;
typedef int64_t                Int64;

/** Unsigned integer type definitions. */
typedef unsigned int           UInt;
typedef unsigned long int      ULInt;
typedef unsigned long long int ULLInt;
typedef uint8_t                UInt8;
typedef uint16_t               UInt16;
typedef uint32_t               UInt32;
typedef uint64_t               UInt64;

/** Integral limits. */
template<typename T = int> requires std::integral<T> constexpr T MinInt(std::numeric_limits<T>::min());
template<typename T = int> requires std::integral<T> constexpr T MaxInt(std::numeric_limits<T>::max());

/***************************************************************************************************************************************************************
* Floating-point Types
***************************************************************************************************************************************************************/

/** Machine epsion and definition of floating point accuracy. */
#define DOUBLE_PRECISION
//#define LONG_DOUBLE_PRECISION

#if defined(DOUBLE_PRECISION)
  typedef double Real;
#elif defined(LONG_DOUBLE_PRECISION)
  typedef long double Real;
#else
  typedef double Real;
  #error No numerical precision defined.
#endif

/** Floating-point limits. */
template<typename T = Real> requires std::floating_point<T> constexpr T Epsilon     (std::numeric_limits<T>::epsilon());
template<typename T = Real> requires std::floating_point<T> constexpr T MinFloat    (std::numeric_limits<T>::min());
template<typename T = Real> requires std::floating_point<T> constexpr T MaxFloat    (std::numeric_limits<T>::max());
template<typename T = Real> requires std::floating_point<T> constexpr T LowestFloat (std::numeric_limits<T>::lowest());
template<typename T = Real> requires std::floating_point<T> constexpr T InfFloat    (std::numeric_limits<T>::infinity());
template<typename T = Real> requires std::floating_point<T> constexpr T QuietNaN    (std::numeric_limits<T>::quiet_NaN());
template<typename T = Real> requires std::floating_point<T> constexpr T SignalNaN   (std::numeric_limits<T>::signaling_NaN());

/** Check if a value is NaN. */
template<typename T>
constexpr bool
isNaN(const T value = T()) { return std::isnan(value); }

/** Check if a value is infinity. */
template<typename T>
constexpr bool
isInfinity(const T value = T()) { return std::isinf(value); }

/***************************************************************************************************************************************************************
* Other Types
***************************************************************************************************************************************************************/

/** Smart pointers/reference wrappers. */
template<class T> using UPtr = std::unique_ptr<T>;
template<class T> using SPtr = std::shared_ptr<T>;
template<class T> using WPtr = std::weak_ptr<T>;
template<class T> using RWpr = std::reference_wrapper<T>;

template<class T> using Option = std::optional<T>;
template<class T1, class T2 = T1> using Pair = std::pair<T1, T2>;

/***************************************************************************************************************************************************************
* Type Checking
***************************************************************************************************************************************************************/

/** Check if two data types are the same. */
template<class T1, class T2>
constexpr bool
isTypeSame() { return std::is_same_v<T1, T2>; }

/** Check if the data types of a sequence of values are the same. */
template<typename T, typename... Ts>
constexpr bool
isTypeHomogeneous() { return (isTypeSame<T, Ts>() && ...); }

/** Check if the data type is a boolean type. */
template<typename T>
constexpr bool
isBoolean(const T = T{}) { return isTypeSame<T, bool>() || isTypeSame<T, Bool>(); }

/** Check if the data type is an integer type. Note: does not include booleans. */
template<typename T>
constexpr bool
isIntegral(const T = T{}) { return std::is_integral_v<T> && !isBoolean<T>(); }

/** Check if the data type is a floating-point type. */
template<typename T>
constexpr bool
isFloatingPoint(const T = T{}) { return std::is_floating_point_v<T>; }

/** Check if the data type is a number type (floating-point or integer type). */
template<typename T>
constexpr bool
isArithmetic(const T = T{}) { return isIntegral<T>() || isFloatingPoint<T>(); }

/** Check if the data type is a number type (floating-point or integer type). */
template<typename T>
constexpr bool
isStringLiteral(const T& = T{}) { return isTypeSame<T, char*>() || isTypeSame<T, std::string>(); }

/** Check if the data type is an enum type. */
template<typename T>
constexpr bool
isEnum(const T = T{}) { return std::is_enum_v<T>; }

/** Check if two enums are the same. */
template<auto T1, auto T2>
requires std::is_enum_v<decltype(T1)> && std::is_enum_v<decltype(T2)>
constexpr bool
isEnumSame() { return T1 == T2; }

/** Get the type category of the given data type (integer, floating-point, etc.). */
template<typename T>
constexpr TypeCategory
GetTypeCategory(const T& = T{})
{
  return isBoolean<T>()       ? TypeCategory::Boolean :
         isIntegral<T>()      ? TypeCategory::Integral :
         isFloatingPoint<T>() ? TypeCategory::FloatingPoint :
         isStringLiteral<T>() ? TypeCategory::String :
                                TypeCategory::Other;
}

/** Get the type of a given iterator. */
template<class It>
using IterType = typename std::iterator_traits<It>::value_type;

/***************************************************************************************************************************************************************
* Type Traits
***************************************************************************************************************************************************************/

/** Check if the data type is an lvalue reference. */
template<typename T>
constexpr bool
isLValueReference() { return std::is_lvalue_reference_v<T>; }

/** Check if the data type is an rvalue reference. */
template<typename T>
constexpr bool
isRValueReference() { return std::is_rvalue_reference_v<T>; }

/** Reference, const, volatile removal */
template<typename T> using RemoveReference = std::remove_reference_t<T>;
template<typename T> using RemoveVolatile  = std::remove_volatile_t<T>;
template<typename T> using RemoveConst     = std::remove_const_t<T>;
template<typename T> using RemoveConstRef  = RemoveConst<RemoveReference<T>>;

/***************************************************************************************************************************************************************
* Static/Dynamic Type Initialisation
***************************************************************************************************************************************************************/

/** Get the initial value for each type category. */
template<typename T, TypeCategory cat = GetTypeCategory<T>()>
constexpr T
StaticInitValue()
{
  return cat == TypeCategory::Boolean       ? static_cast<T>(false) :
         cat == TypeCategory::Integral      ? static_cast<T>(-1) :
         cat == TypeCategory::FloatingPoint ? static_cast<T>(0.0) :
         cat == TypeCategory::Other         ? T{} :
         throw std::invalid_argument("The passed type does not qualify for static initialisation.");
}

/** Get the initial value for each type category. */
template<typename T, TypeCategory cat = GetTypeCategory<T>()>
T
DynamicInitValue()
{
  return cat == TypeCategory::Boolean       ? static_cast<T>(false) :
         cat == TypeCategory::Integral      ? static_cast<T>(-1) :
         cat == TypeCategory::FloatingPoint ? static_cast<T>(0.0) :
         cat == TypeCategory::String        ? static_cast<T>('\0') :
         cat == TypeCategory::Other         ? T{} :
         throw std::invalid_argument("The passed type does not qualify for dynamic initialisation.");
}

/***************************************************************************************************************************************************************
* Type Casting
***************************************************************************************************************************************************************/

/** Convert an enum to its corresponding integer value. */
template <typename T>
requires std::is_enum_v<T>
constexpr auto EnumToInt(T e) noexcept { return static_cast<std::underlying_type_t<T>>(e); }

/** Convert an rvalue reference to an lvalue reference. */
template <typename T>
constexpr T& Unmove(T&& t) noexcept { return t; }

}
