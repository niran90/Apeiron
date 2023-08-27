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
#include "../include/Global.h"
#include "ApeironTest.h"

#ifdef DEBUG_MODE

namespace aprn {

/***************************************************************************************************************************************************************
* Sequences and Series
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, GetInteger)
{
  constexpr size_t size = 100;
  FOR(i, size) EXPECT_EQ(GetInteger(std::make_integer_sequence<int, size>{}, i), i);
}

/***************************************************************************************************************************************************************
* General Arithmetic Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, Sum)
{
  constexpr int intSum = Sum(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
  EXPECT_EQ(intSum, 55);

  constexpr Real floatSum = Sum(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0);
  EXPECT_DOUBLE_EQ(floatSum, 55.0);
}

TEST_F(ApeironTest, Product)
{
  constexpr int intProd = Product(1, 2, 3, 4, 5);
  EXPECT_EQ(intProd, 120);

  constexpr Real floatProd = Product(1.0, 2.0, 3.0, 4.0, 5.0);
  EXPECT_DOUBLE_EQ(floatProd, 120.0);
}

TEST_F(ApeironTest, Divide)
{
  EXPECT_EQ(Divide(20, 5), Four);
  EXPECT_THROW(Divide(1, 0), std::logic_error);

  EXPECT_EQ(Divide(20.0, 5.0), Four);
  EXPECT_THROW(Divide(1.0, 0.0), std::logic_error);
}

TEST_F(ApeironTest, Modulo)
{
  EXPECT_EQ(Modulo(20, 5), 0);
  EXPECT_EQ(Modulo(20, 3), 2);
  EXPECT_THROW(Modulo(1, 0), std::logic_error);

  EXPECT_DOUBLE_EQ(Modulo(20.0, 5.0), Zero);
  EXPECT_DOUBLE_EQ(Modulo(20.0, 3.0), Two);
  EXPECT_DOUBLE_EQ(Modulo(20.1234, 3.0), 2.1234);
  EXPECT_THROW(Modulo(1.0, 0.0), std::logic_error);
}

/***************************************************************************************************************************************************************
* Combinatorial Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, Factorial)
{
  EXPECT_EQ(Factorial(0), 1);
  EXPECT_EQ(Factorial(1), 1);
  EXPECT_EQ(Factorial(5), 120);
}

TEST_F(ApeironTest, FactorialQuotient)
{
  EXPECT_EQ(FactorialQuotient(1, 0), 1);
  EXPECT_EQ(FactorialQuotient(2, 1), 2);
  EXPECT_EQ(FactorialQuotient(2, 2), 1);
  EXPECT_EQ(FactorialQuotient(5, 3), 20);
  EXPECT_THROW(FactorialQuotient(2, 3), std::logic_error);
}

TEST_F(ApeironTest, Choose)
{
  EXPECT_EQ(Choose(3, 0), 1);
  EXPECT_EQ(Choose(3, 1), 3);
  EXPECT_EQ(Choose(3, 2), 3);
  EXPECT_EQ(Choose(3, 3), 1);
  EXPECT_THROW(Choose(3, 4), std::logic_error);
}

/***************************************************************************************************************************************************************
* Power Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, iPow)
{
  EXPECT_EQ(iPow(2.0, 1), 2.0);
  EXPECT_EQ(iPow(2.0, 2), 4.0);
  EXPECT_EQ(iPow(2.0, 3), 8.0);
  EXPECT_THROW(iPow(2.0, 31), std::logic_error);
}

}

#endif