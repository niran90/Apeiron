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

#include <gtest/gtest.h>
#include "../include/NumericContainer.h"
#include "../../DataContainer/include/Array.h"

#ifdef DEBUG_MODE

constexpr size_t ContainerSize = 50;

namespace aprn {
namespace Test {

/***************************************************************************************************************************************************************
* NumericContainer Test Classes
***************************************************************************************************************************************************************/
template<typename T, class D>
class NumericContainer : public detail::NumericContainer<T, D>
{
protected:
  constexpr NumericContainer() = default;

public:
  constexpr D& Derived() noexcept { return static_cast<D&>(*this); }
  constexpr const D& Derived() const noexcept { return static_cast<const D&>(*this); }
};

template<typename T, size_t N>
struct StaticNumericContainer : public StaticArray<T, N>,
                               public NumericContainer<T, StaticNumericContainer<T, N>>
{

};

template<typename T>
struct DynamicNumericContainer : public DynamicArray<T>,
                                public NumericContainer<T, DynamicNumericContainer<T>>
{
  DynamicNumericContainer(const size_t size) : DynamicArray<T>::DynamicArray(size) {}
};

/***************************************************************************************************************************************************************
* NumericContainer Test Fixture
***************************************************************************************************************************************************************/
class NumericContainerTest : public testing::Test
{
public:
  Random<int>  RandomInt;
  Random<Real> RandomReal;

  StaticNumericContainer<int, ContainerSize>   IntStaticContainer;
  StaticNumericContainer <Real, ContainerSize> RealStaticContainer;

  DynamicNumericContainer<int>   IntDynamicContainer;
  DynamicNumericContainer <Real> RealDynamicContainer;

  NumericContainerTest()
    : RandomInt(-10, 10), RandomReal(-Ten, Ten), IntDynamicContainer(ContainerSize), RealDynamicContainer(ContainerSize) {}

  void SetUp() override
  {
    FOR_EACH(entry, IntStaticContainer)    entry = RandomInt();
    FOR_EACH(entry, IntDynamicContainer)   entry = RandomInt();
    FOR_EACH(entry, RealStaticContainer)   entry = RandomReal();
    FOR_EACH(entry, RealDynamicContainer) entry = RandomReal();
  }
};

}
}

#endif