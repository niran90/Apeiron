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
#include "../../../include/Random.h"
#include "../include/Array.h"

#ifdef DEBUG_MODE

constexpr size_t ContainerSize = 50;

namespace aprn
{

class ArrayTest : public testing::Test
{
public:
  StaticArray<bool , ContainerSize> BoolStaticArray;
  StaticArray<int  , ContainerSize> IntStaticArray;
  StaticArray<Real, ContainerSize>  RealStaticArray;

  DynamicArray<Bool> BoolDynamicArray;
  DynamicArray<int>  IntDynamicArray;
  DynamicArray<Real> RealDynamicArray;

  ArrayTest()
    : BoolDynamicArray(ContainerSize), IntDynamicArray(ContainerSize), RealDynamicArray(ContainerSize) {}

  void SetUp() override
  {

  }
};

TEST_F(ArrayTest, Initialisation)
{
  // Initilise constexpr StaticArray objects.
  constexpr StaticArray<bool, ContainerSize> ConstexprBoolArray;
  constexpr StaticArray<int, ContainerSize>  ConstexprIntArray(2);
  constexpr StaticArray<Real, 5> ConstexprRealArray1{One, Two, Three, Four, Five};
  constexpr StaticArray<Real, 3> ConstexprRealArray2(ConstexprRealArray1.begin(), ConstexprRealArray1.begin() + 3);

  FOR(i, 5)
  {
    EXPECT_EQ(ConstexprRealArray1[i], Real(i + 1));
    if(i < 3) { EXPECT_EQ(ConstexprRealArray2[i], Real(i + 1)); }
  }

  // Test default initialisation of arrays.
  FOR(i, ContainerSize)
  {
    // Static arrays
    EXPECT_EQ(BoolStaticArray[i], false);
    EXPECT_EQ(IntStaticArray[i], -1);
    EXPECT_EQ(RealStaticArray[i], Zero);

    // Dynamic arrays
    EXPECT_EQ(BoolDynamicArray[i], false);
    EXPECT_EQ(IntDynamicArray[i], -1);
    EXPECT_EQ(RealDynamicArray[i], Zero);
  }
}

TEST_F(ArrayTest, SubscriptOperator)
{
  FOR(i, ContainerSize)
  {
    // Static arrays
    BoolStaticArray[i] = true;
    IntStaticArray[i] = 1;
     RealStaticArray[i] = Zero;
    EXPECT_EQ(BoolStaticArray[i], true);
    EXPECT_EQ(IntStaticArray[i], 1);
    EXPECT_EQ(RealStaticArray[i], Zero);

    // Dynamic arrays
    BoolDynamicArray[i] = true;
    IntDynamicArray[i] = 1;
     RealDynamicArray[i] = Zero;
    EXPECT_EQ(BoolDynamicArray[i], true);
    EXPECT_EQ(IntDynamicArray[i], 1);
    EXPECT_EQ(RealDynamicArray[i], Zero);
  }
}

TEST_F(ArrayTest, AssignmentOperator)
{
  const bool rand_bool = Random<bool>()();
  const int rand_int = Random(-100, 100)();
  const Real rand_float = Random(-100.0, 100.0)();

  BoolStaticArray = rand_bool;
  IntStaticArray = rand_int;
   RealStaticArray = rand_float;

  BoolDynamicArray = rand_bool;
  IntDynamicArray = rand_int;
   RealDynamicArray = rand_float;

  FOR(i, ContainerSize)
  {
    // Static arrays
    EXPECT_EQ(BoolStaticArray[i], rand_bool);
    EXPECT_EQ(IntStaticArray[i], rand_int);
    EXPECT_EQ(RealStaticArray[i], rand_float);

    // Dynamic arrays
    EXPECT_EQ(BoolDynamicArray[i], rand_bool);
    EXPECT_EQ(IntDynamicArray[i], rand_int);
    EXPECT_EQ(RealDynamicArray[i], rand_float);
  }
}

}

#endif

