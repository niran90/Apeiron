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

#include <algorithm>
#include <tuple>

namespace aprn {

/***************************************************************************************************************************************************************
* Conditional Helper Functions
***************************************************************************************************************************************************************/
namespace detail {

/** OneOf helper struct for a variable number of arguments. */
template<typename ...Ts>
struct VariadicOneOfHelper
{
   constexpr VariadicOneOfHelper(Ts... values) : Values(std::move(values)...) {}

   template<typename T>
   [[nodiscard]] friend constexpr bool
   operator==(T lhs, const VariadicOneOfHelper& rhs) noexcept { return std::apply([&](auto... vals){ return ((lhs == vals) || ...); }, rhs.Values); }

   template<typename T>
   [[nodiscard]] friend constexpr bool
   operator!=(T lhs, const VariadicOneOfHelper& rhs) noexcept { return !operator==(lhs, rhs); }

   std::tuple<Ts...> Values;
};

/** OneOf helper struct for entries of an STL container type. */
template<typename Cont>
struct ContainerOneOfHelper
{
   constexpr ContainerOneOfHelper(Cont const& values) : Values(values) {}

   template<typename T>
   [[nodiscard]] friend constexpr bool
   operator==(T&& lhs, ContainerOneOfHelper&& rhs) noexcept { return std::any_of(cbegin(rhs.Values), cend(rhs.Values), [&](auto val){ return lhs == val; }); }

   template<typename T>
   [[nodiscard]] friend constexpr bool
   operator!=(T&& lhs, ContainerOneOfHelper&& rhs) noexcept { return !operator==(std::move(lhs), std::move(rhs)); }

   const Cont& Values;
};

template<typename T>
auto isContainer(int) -> decltype(cbegin(std::declval<T>()) == cend(std::declval<T>()), std::true_type{});

template<typename T>
std::false_type isContainer(...);

template<typename T>
constexpr bool isContainer_v = decltype(isContainer<T>(0))::value;

}//detail

/***************************************************************************************************************************************************************
* Conditionals
***************************************************************************************************************************************************************/

/** Check for equality with 'any of' the entries in a list of arguments or entries in an STL container. */
template<typename ...Ts>
[[nodiscard]] constexpr auto
OneOf(Ts&&... values)
{
   using namespace detail;
   if constexpr(sizeof...(Ts) == 1 && isContainer_v<std::tuple_element_t<0, std::tuple<Ts...>>>) return ContainerOneOfHelper(std::forward<Ts>(values)...);
   else return VariadicOneOfHelper(std::forward<Ts>(values)...);
}

}
