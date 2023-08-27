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
#include "../include/Random.h"
#include "ApeironTest.h"

#ifdef DEBUG_MODE

namespace aprn{

TEST_F(ApeironTest, isEqual)
{
  const Real rela_tol = RelativeTolerance;
  const Real zero_tol = ZeroTolerance;

  // Zero toleranced equality (Reals only)
  EXPECT_TRUE(isEqual(Zero, Zero));
  EXPECT_TRUE(isEqual(Zero, 0.999 * zero_tol));
  EXPECT_TRUE(isEqual(0.999 * zero_tol, Zero));
  EXPECT_FALSE(isEqual(Zero, zero_tol));
  EXPECT_FALSE(isEqual(zero_tol, Zero));
  EXPECT_FALSE(isEqual(Zero, One));

  FOR(i, 1000)
  {
    const int random_int = RandomInt();
    const Real random_float = RandomReal();

    // Exact precision equality
    EXPECT_TRUE(isEqual<true>(random_int, random_int));
    EXPECT_TRUE(isEqual<true>(random_float, random_float));
    EXPECT_FALSE(isEqual<true>(random_int, random_int + 1));
    EXPECT_FALSE(isEqual<true>(random_float, random_float + Small));

    // Relative toleranced equality (Reals only)
    EXPECT_TRUE(isEqual(random_float, random_float));
    EXPECT_TRUE(isEqual(random_float, random_float + 0.9 * rela_tol * random_float));
    EXPECT_TRUE(isEqual(random_float + 0.9 * rela_tol * random_float, random_float));
    EXPECT_TRUE(isEqual(random_float, random_float - 0.9 * rela_tol * random_float));
    EXPECT_TRUE(isEqual(random_float - 0.9 * rela_tol * random_float, random_float));
    EXPECT_FALSE(isEqual(random_float, random_float + 1.1 * rela_tol * random_float));
    EXPECT_FALSE(isEqual(random_float + 1.1 * rela_tol * random_float, random_float));
    EXPECT_FALSE(isEqual(random_float, random_float - 1.1 * rela_tol * random_float));
    EXPECT_FALSE(isEqual(random_float - 1.1 * rela_tol * random_float, random_float));
    EXPECT_FALSE(isEqual(random_float, random_float + One));
  }
}

TEST_F(ApeironTest, isLess)
{
  const Real rela_tol = RelativeTolerance;
  const Real zero_tol = ZeroTolerance;

  // Zero toleranced lesser than (Reals only)
  EXPECT_TRUE(isLess(Zero, One));
  EXPECT_TRUE(isLess(Zero, zero_tol));
  EXPECT_FALSE(isLess(Zero, 0.999 * zero_tol));
  EXPECT_FALSE(isLess(Zero, Zero));

  FOR(i, 1000)
  {
    const int random_int = RandomInt();
    const Real random_float = RandomReal();

    // Exact precision lesser than
    EXPECT_TRUE(isLess<true>(random_int, random_int + 1));
    EXPECT_TRUE(isLess<true>(random_float, random_float + Small));
    EXPECT_FALSE(isLess<true>(random_int, random_int));
    EXPECT_FALSE(isLess<true>(random_float, random_float));

    // Relative toleranced lesser than (Reals only)
    EXPECT_TRUE(isLess(random_float, random_float + 1.1 * rela_tol * Abs(random_float)));
    EXPECT_TRUE(isLess(random_float - 1.1 * rela_tol * Abs(random_float), random_float));
    EXPECT_TRUE(isLess(random_float, random_float + One));
    EXPECT_FALSE(isLess(random_float, random_float));
    EXPECT_FALSE(isLess(random_float, random_float + 0.9 * rela_tol * Abs(random_float)));
    EXPECT_FALSE(isLess(random_float - 0.9 * rela_tol * Abs(random_float), random_float));
  }
}

TEST_F(ApeironTest, isLessEqual)
{
  const Real rela_tol = RelativeTolerance;
  const Real zero_tol = ZeroTolerance;

  // Zero toleranced lesser than or equal to (Reals only)
  EXPECT_TRUE(isLessEqual(Zero, Zero));
  EXPECT_TRUE(isLessEqual(Zero, One));
  EXPECT_TRUE(isLessEqual(Zero, 0.999 * zero_tol));
  EXPECT_FALSE(isLessEqual(1.111 * zero_tol, Zero));

  FOR(i, 1000)
  {
    const int  random_int = RandomInt();
    const Real random_real = RandomReal();

    // Exact precision lesser than or equal to
    EXPECT_TRUE(isLessEqual<true>(random_int, random_int));
    EXPECT_TRUE(isLessEqual<true>(random_real, random_real));
    EXPECT_TRUE(isLessEqual<true>(random_int, random_int + 1));
    EXPECT_TRUE(isLessEqual<true>(random_real, random_real + Small));

    // Relative toleranced lesser than or equal to (Floats only)
    EXPECT_TRUE(isLessEqual(random_real, random_real));
    EXPECT_TRUE(isLessEqual(random_real, random_real + One));
    EXPECT_TRUE(isLessEqual(random_real, random_real + 0.9 * rela_tol * Abs(random_real)));
    EXPECT_TRUE(isLessEqual(random_real - 1.1 * rela_tol * Abs(random_real), random_real));
    EXPECT_FALSE(isLessEqual(random_real + 1.1 * rela_tol * Abs(random_real), random_real));
  }
}

TEST_F(ApeironTest, isLarger)
{
  FOR(i, 1000)
  {
    const int random_int = RandomInt();
    const Real random_float = RandomReal();

    EXPECT_EQ(isLarger<true>(random_int, random_int), !isLessEqual<true>(random_int, random_int));
    EXPECT_EQ(isLarger<true>(random_float, random_float), !isLessEqual<true>(random_float, random_float));
    EXPECT_EQ(isLarger(random_float, random_float), !isLessEqual(random_float, random_float));
  }
}

TEST_F(ApeironTest, isLargerEqual)
{
  FOR(i, 1000)
  {
    const int random_int = RandomInt();
    const Real random_float = RandomReal();

    EXPECT_EQ(isLargerEqual<true>(random_int, random_int), !isLess<true>(random_int, random_int));
    EXPECT_EQ(isLargerEqual<true>(random_float, random_float), !isLess<true>(random_float, random_float));
    EXPECT_EQ(isLargerEqual(random_float, random_float), !isLess(random_float, random_float));
  }
}

TEST_F(ApeironTest, isBounded)
{
  FOR(i, 1000)
  {
    const int min_int = -Abs(RandomInt());
    const int max_int = min_int + 1 + Abs(RandomInt());
    const Real min_float = -Abs(RandomReal());
    const Real max_float = min_float + One + Abs(RandomReal());

    // Test left bounds
    EXPECT_TRUE(isBounded(min_int, min_int, max_int));
    EXPECT_TRUE(isBounded<true>(min_int, min_int, max_int));
    EXPECT_FALSE(isBounded<false>(min_int, min_int, max_int));
    EXPECT_TRUE(isBounded(min_float, min_float, max_float));
    EXPECT_TRUE(isBounded<true>(min_float, min_float, max_float));
    EXPECT_FALSE(isBounded<false>(min_float, min_float, max_float));

    // Test right bounds
    EXPECT_FALSE(isBounded(max_int, min_int, max_int));
    EXPECT_FALSE((isBounded<true, false>(max_int, min_int, max_int)));
    EXPECT_TRUE((isBounded<true, true>(max_int, min_int, max_int)));
    EXPECT_FALSE(isBounded(max_float, min_float, max_float));
    EXPECT_FALSE((isBounded<true, false>(max_float, min_float, max_float)));
    EXPECT_TRUE((isBounded<true, true>(max_float, min_float, max_float)));

    // Test inclusion/exclusion
//    const int random_int = RandomInt();
//    const Real random_float = RandomFloat();
//
//    EXPECT_EQ((isBounded<false, false, true>(random_int, min_int, max_int)), ());
//    EXPECT_EQ(isBounded<true>(random_int, min_int, max_int));
//    EXPECT_EQ(isBounded<false>(random_int, min_int, max_int));
//    EXPECT_EQ(isBounded(random_float, min_float, max_float));
//    EXPECT_EQ(isBounded<true>(random_float, min_float, max_float));
//    EXPECT_EQ(isBounded<false>(random_float, min_float, max_float));

  }
}

}

#endif