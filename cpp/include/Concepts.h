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

#include "Types.h"

#include <type_traits>
#include <concepts>

namespace aprn {

/***************************************************************************************************************************************************************
* Type Concepts
***************************************************************************************************************************************************************/
template<typename T> concept Enumerator    = isEnum<T>();
template<typename T> concept Boolean       = isBoolean<T>();
template<typename T> concept Integral      = isIntegral<T>();
template<typename T> concept FloatingPoint = isFloatingPoint<T>();
template<typename T> concept Arithmetic    = isArithmetic<T>();
template<typename T> concept StringLiteral = isStringLiteral<T>();

/***************************************************************************************************************************************************************
* Value Category Concepts
***************************************************************************************************************************************************************/
template<typename T> concept LValue        = isLValueReference<T>();
template<typename T> concept RValue        = isRValueReference<T>();

/***************************************************************************************************************************************************************
* Iterator Concepts
***************************************************************************************************************************************************************/
template<typename T> concept CharIterator  = isTypeSame<typename std::iterator_traits<T>::value_type, char>();

}
