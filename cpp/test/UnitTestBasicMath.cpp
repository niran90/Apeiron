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

/***************************************************************************************************************************************************************
* Counting Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, nDigits)
{
   const auto max_number = 1e5;
   Random<size_t> RandomNumber(0, max_number);

   FOR(i, max_number)
   {
      const auto num = RandomNumber();

      // Manually count digits.
      size_t n_digits(0);
      if     (num <= 9)                     n_digits = 1;
      else if(10 <= num && num <= 99)       n_digits = 2;
      else if(100 <= num && num <= 999)     n_digits = 3;
      else if(1000 <= num && num <= 9999)   n_digits = 4;
      else if(10000 <= num && num <= 99999) n_digits = 5;
      else if(100000 == num)                n_digits = 6;
      else EXIT("Ensure that the maximum number of digits for this test is 6.")

      EXPECT_EQ(nDigits(num), n_digits);
   }
}

/***************************************************************************************************************************************************************
* Min, Max, and Bounding Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, Min)
{
  const int integer = RandomInt();
  EXPECT_EQ(Min(integer, integer + 1), integer);
  EXPECT_EQ(Min(integer, integer - 1), integer - 1);

  const Real floating = RandomReal();
  EXPECT_EQ(Min(floating, floating + One), floating);
  EXPECT_EQ(Min(floating, floating - One), floating - One);
}

TEST_F(ApeironTest, Max)
{
  const int integer = RandomInt();
  EXPECT_EQ(Max(integer, integer + 1), integer + 1);
  EXPECT_EQ(Max(integer, integer - 1), integer);

  const Real floating = RandomReal();
  EXPECT_EQ(Max(floating, floating + One), floating + One);
  EXPECT_EQ(Max(floating, floating - One), floating);
}

TEST_F(ApeironTest, Clipped)
{
  const int minInt = -Abs(RandomInt());
  const int maxInt = Abs(RandomInt()) + 1;
  const int integer = RandomInt();
  EXPECT_GE(Clipped(integer, minInt, maxInt), minInt);
  EXPECT_LE(Clipped(integer, minInt, maxInt), maxInt);

  const Real minFloat = -Abs(RandomReal());
  const Real maxFloat = Abs(RandomReal()) + One;
  const Real floating = RandomReal();
  EXPECT_GE(Clipped(floating, minFloat, maxFloat), minFloat);
  EXPECT_LE(Clipped(floating, minFloat, maxFloat), maxFloat);
}

TEST_F(ApeironTest, MinEntry)
{
  int minInt = MaxInt<>;
  FOR_EACH(entry, IntArray) minInt = Min(minInt, entry);
  EXPECT_EQ(MinEntry(IntArray.begin(), IntArray.end()), minInt);

  Real minFloat = MaxFloat<>;
  FOR_EACH(entry, RealArray) minFloat = Min(minFloat, entry);
  EXPECT_EQ(MinEntry(RealArray.begin(), RealArray.end()), minFloat);
}

TEST_F(ApeironTest, MaxEntry)
{
  int maxInt = MinInt<>;
  FOR_EACH(entry, IntArray) maxInt = Max(maxInt, entry);
  EXPECT_EQ(MaxEntry(IntArray.begin(), IntArray.end()), maxInt);

  Real maxFloat = MinFloat<>;
  FOR_EACH(entry, RealArray) maxFloat = Max(maxFloat, entry);
  EXPECT_EQ(MaxEntry(RealArray.begin(), RealArray.end()), maxFloat);
}

TEST_F(ApeironTest, MinMaxEntries)
{
  int minInt = MaxInt<>;
  int maxInt = MinInt<>;
  FOR_EACH(entry, IntArray) minInt = Min(minInt, entry);
  FOR_EACH(entry, IntArray) maxInt = Max(maxInt, entry);
  EXPECT_EQ(MinMaxEntries(IntArray.begin(), IntArray.end()).first, minInt);
  EXPECT_EQ(MinMaxEntries(IntArray.begin(), IntArray.end()).second, maxInt);

  Real minFloat = MaxFloat<>;
  Real maxFloat = MinFloat<>;
  FOR_EACH(entry, RealArray) minFloat = Min(minFloat, entry);
  FOR_EACH(entry, RealArray) maxFloat = Max(maxFloat, entry);
  EXPECT_EQ(MinMaxEntries(RealArray.begin(), RealArray.end()).first, minFloat);
  EXPECT_EQ(MinMaxEntries(RealArray.begin(), RealArray.end()).second, maxFloat);
}

TEST_F(ApeironTest, BoundEntries)
{
  const int minInt = -Abs(RandomInt());
  const int maxInt = Abs(RandomInt()) + 1;
  BoundEntries(IntArray.begin(), IntArray.end(), minInt, maxInt);
  FOR_EACH(entry, IntArray)
  {
    EXPECT_GE(entry, minInt);
    EXPECT_LE(entry, maxInt);
  }

  const Real minFloat = -Abs(RandomReal());
  const Real maxFloat = Abs(RandomReal()) + One;
  BoundEntries(RealArray.begin(), RealArray.end(), minFloat, maxFloat);
  FOR_EACH(entry, RealArray)
  {
    EXPECT_GE(entry, minFloat);
    EXPECT_LE(entry, maxFloat);
  }
}

/***************************************************************************************************************************************************************
* Signum and Absolute Value Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, Sgn)
{
  const int negativeInt = -Abs(RandomInt()) - 1;
  const int positiveInt = Abs(RandomInt()) + 1;
  EXPECT_EQ(Sgn(negativeInt), -1);
  EXPECT_EQ(Sgn(positiveInt), 1);
  EXPECT_EQ(Sgn(0), 1);
  EXPECT_EQ(Sgn(0, 1), 1);
  EXPECT_EQ(Sgn(0, 0), 0);
  EXPECT_EQ(Sgn(0, -1), -1);

  const Real negativeFloat = -Abs(RandomReal()) - One;
  const Real positiveFloat = Abs(RandomReal()) + One;
  EXPECT_EQ(Sgn(negativeFloat), -One);
  EXPECT_EQ(Sgn(positiveFloat), One);
  EXPECT_EQ(Sgn(Zero), One);
  EXPECT_EQ(Sgn(Zero, 1), One);
  EXPECT_EQ(Sgn(Zero, 0), Zero);
  EXPECT_EQ(Sgn(Zero, -1), -One);
}

TEST_F(ApeironTest, Abs)
{
  const int negativeInt = -Abs(RandomInt()) - 1;
  const int positiveInt = Abs(RandomInt()) + 1;
  EXPECT_EQ(Abs(negativeInt), -negativeInt);
  EXPECT_EQ(Abs(positiveInt), positiveInt);

  const Real negativeFloat = -Abs(RandomReal()) - One;
  const Real positiveFloat = Abs(RandomReal()) + One;
  EXPECT_EQ(Abs(negativeFloat), -negativeFloat);
  EXPECT_EQ(Abs(positiveFloat), positiveFloat);
}

/***************************************************************************************************************************************************************
* Numerical Rounding Functions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, Floor)
{
  EXPECT_EQ(Floor(-1.234), -Two);
  EXPECT_EQ(Floor(-1.678), -Two);
  EXPECT_EQ(Floor(-Two), -Two);

  EXPECT_EQ(Floor(1.234), One);
  EXPECT_EQ(Floor(1.999), One);
  EXPECT_EQ(Floor(One), One);
}

TEST_F(ApeironTest, Ceil)
{
  EXPECT_EQ(Ceil(-1.234), -One);
  EXPECT_EQ(Ceil(-1.999), -One);
  EXPECT_EQ(Ceil(-Two), -Two);

  EXPECT_EQ(Ceil(1.234), Two);
  EXPECT_EQ(Ceil(1.999), Two);
  EXPECT_EQ(Ceil(One), One);
}

/***************************************************************************************************************************************************************
* Mathematical Conversions
***************************************************************************************************************************************************************/
TEST_F(ApeironTest, DegToRad)
{
  EXPECT_DOUBLE_EQ(DegToRad(30.0), SixthPi);
  EXPECT_DOUBLE_EQ(DegToRad(60.0), ThirdPi);
  EXPECT_DOUBLE_EQ(DegToRad(90.0), HalfPi);
  EXPECT_DOUBLE_EQ(DegToRad(180.0), Pi);
}

TEST_F(ApeironTest, RadToDeg)
{
  EXPECT_DOUBLE_EQ(RadToDeg(SixthPi), 30.0);
  EXPECT_DOUBLE_EQ(RadToDeg(ThirdPi), 60.0);
  EXPECT_DOUBLE_EQ(RadToDeg(HalfPi), 90.0);
  EXPECT_DOUBLE_EQ(RadToDeg(Pi), 180.0);
}

}

#endif