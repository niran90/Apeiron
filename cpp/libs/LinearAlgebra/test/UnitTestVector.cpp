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

#include "../../../include/Global.h"
#include "../include/Vector.h"
#include "../include/VectorOperations.h"

#ifdef DEBUG_MODE

constexpr size_t ContainerSize = 50;

namespace aprn {

/***************************************************************************************************************************************************************
* Vector Test Fixture
***************************************************************************************************************************************************************/
class VectorTest : public testing::Test
{
 public:
   Random<int>  RandomInt;
   Random<Real> RandomReal;

   StaticVector<int, ContainerSize>  IntStaticVector;
   StaticVector<Real, ContainerSize> RealStaticVector;
   DynamicVector<int>  IntDynamicVector;
   DynamicVector<Real> RealDynamicVector;

   StaticVector<int, ContainerSize>  IntStaticVectorTest;
   StaticVector<Real, ContainerSize> RealStaticVectorTest;
   DynamicVector<int>  IntDynamicVectorTest;
   DynamicVector<Real> RealDynamicVectorTest;

   VectorTest()
      : RandomInt(-10, 10), RandomReal(-Ten, Ten), IntDynamicVector(ContainerSize), RealDynamicVector(ContainerSize),
        IntDynamicVectorTest(ContainerSize), RealDynamicVectorTest(ContainerSize) {}

   void SetUp() override
   {
      IntStaticVector.ResetRandomiser(-10, 10);
      IntDynamicVector.ResetRandomiser(-10, 10);
      RealStaticVector.ResetRandomiser(-Ten, Ten);
      RealDynamicVector.ResetRandomiser(-Ten, Ten);

      IntStaticVector.Randomise();
      IntDynamicVector.Randomise();
      RealStaticVector.Randomise();
      RealDynamicVector.Randomise();

      IntStaticVectorTest = IntStaticVector;
      IntDynamicVectorTest = IntDynamicVector;
      RealStaticVectorTest = RealStaticVector;
      RealDynamicVectorTest = RealDynamicVector;
   }
};

/***************************************************************************************************************************************************************
* Test Arithmetic Operations
***************************************************************************************************************************************************************/
TEST_F(VectorTest, ScalarPlus)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  IntStaticVectorTest = IntStaticVector + random_int;
  IntDynamicVectorTest = IntDynamicVector + random_int;
   RealStaticVectorTest = RealStaticVector + random_float;
   RealDynamicVectorTest = RealDynamicVector + random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] + random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] + random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] + random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] + random_float);
  }
}

TEST_F(VectorTest, ScalarMinus)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  IntStaticVectorTest = IntStaticVector - random_int;
  IntDynamicVectorTest = IntDynamicVector - random_int;
   RealStaticVectorTest = RealStaticVector - random_float;
   RealDynamicVectorTest = RealDynamicVector - random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] - random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] - random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] - random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] - random_float);
  }
}

TEST_F(VectorTest, ScalarMultiply)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  // Test pre-multiplication
  IntStaticVectorTest = random_int * IntStaticVector;
  IntDynamicVectorTest = random_int * IntDynamicVector;
   RealStaticVectorTest = random_float * RealStaticVector;
   RealDynamicVectorTest = random_float * RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], random_int * IntStaticVector[i]);
    EXPECT_EQ(IntDynamicVectorTest[i], random_int * IntDynamicVector[i]);
    EXPECT_EQ(RealStaticVectorTest[i], random_float * RealStaticVector[i]);
    EXPECT_EQ(RealDynamicVectorTest[i], random_float * RealDynamicVector[i]);
  }

  // Test post-multiplication
  IntStaticVectorTest = IntStaticVector * random_int;
  IntDynamicVectorTest = IntDynamicVector * random_int;
   RealStaticVectorTest = RealStaticVector * random_float;
   RealDynamicVectorTest = RealDynamicVector * random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], random_int * IntStaticVector[i]);
    EXPECT_EQ(IntDynamicVectorTest[i], random_int * IntDynamicVector[i]);
    EXPECT_EQ(RealStaticVectorTest[i], random_float * RealStaticVector[i]);
    EXPECT_EQ(RealDynamicVectorTest[i], random_float * RealDynamicVector[i]);
  }
}

TEST_F(VectorTest, ScalarDivide)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  if(random_int == 0) random_int++;
  if(isEqual(random_float, Zero)) random_float++;

  IntStaticVectorTest = IntStaticVector / random_int;
  IntDynamicVectorTest = IntDynamicVector / random_int;
   RealStaticVectorTest = RealStaticVector / random_float;
   RealDynamicVectorTest = RealDynamicVector / random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] / random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] / random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] / random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] / random_float);
  }
}

TEST_F(VectorTest, ScalarPlusEqual)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  IntStaticVectorTest += random_int;
  IntDynamicVectorTest += random_int;
   RealStaticVectorTest += random_float;
   RealDynamicVectorTest += random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] + random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] + random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] + random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] + random_float);
  }
}

TEST_F(VectorTest, ScalarMinusEqual)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  IntStaticVectorTest -= random_int;
  IntDynamicVectorTest -= random_int;
   RealStaticVectorTest -= random_float;
   RealDynamicVectorTest -= random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] - random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] - random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] - random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] - random_float);
  }
}

TEST_F(VectorTest, ScalarMultiplyEqual)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  IntStaticVectorTest *= random_int;
  IntDynamicVectorTest *= random_int;
   RealStaticVectorTest *= random_float;
   RealDynamicVectorTest *= random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] * random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] * random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] * random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] * random_float);
  }
}

TEST_F(VectorTest, ScalarDivideEqual)
{
  int random_int = RandomInt();
  Real random_float = RandomReal();

  if(random_int == 0) random_int++;
  if(isEqual(random_float, Zero)) random_float++;

  IntStaticVectorTest /= random_int;
  IntDynamicVectorTest /= random_int;
   RealStaticVectorTest /= random_float;
   RealDynamicVectorTest /= random_float;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], IntStaticVector[i] / random_int);
    EXPECT_EQ(IntDynamicVectorTest[i], IntDynamicVector[i] / random_int);
    EXPECT_EQ(RealStaticVectorTest[i], RealStaticVector[i] / random_float);
    EXPECT_EQ(RealDynamicVectorTest[i], RealDynamicVector[i] / random_float);
  }
}

TEST_F(VectorTest, VectorPlus)
{
  IntStaticVectorTest = IntStaticVector + IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector + IntDynamicVector;
   RealStaticVectorTest = RealStaticVector + RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector + RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 2 * IntStaticVector[i]);
    EXPECT_EQ(IntDynamicVectorTest[i], 2 * IntDynamicVector[i]);
    EXPECT_EQ(RealStaticVectorTest[i], Two * RealStaticVector[i]);
    EXPECT_EQ(RealDynamicVectorTest[i], Two * RealDynamicVector[i]);
  }
}

TEST_F(VectorTest, VectorMinus)
{
  IntStaticVectorTest = IntStaticVector - IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector - IntDynamicVector;
   RealStaticVectorTest = RealStaticVector - RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector - RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 0);
    EXPECT_EQ(IntDynamicVectorTest[i], 0);
    EXPECT_EQ(RealStaticVectorTest[i], Zero);
    EXPECT_EQ(RealDynamicVectorTest[i], Zero);
  }
}

TEST_F(VectorTest, VectorMultiply)
{
  IntStaticVectorTest = IntStaticVector * IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector * IntDynamicVector;
   RealStaticVectorTest = RealStaticVector * RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector * RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], iPow(IntStaticVector[i], 2));
    EXPECT_EQ(IntDynamicVectorTest[i], iPow(IntDynamicVector[i], 2));
    EXPECT_EQ(RealStaticVectorTest[i], iPow(RealStaticVector[i], 2));
    EXPECT_EQ(RealDynamicVectorTest[i], iPow(RealDynamicVector[i], 2));
  }
}

TEST_F(VectorTest, VectorDivide)
{
  FOR_EACH(entry, IntStaticVector) if(entry == 0) entry++;
  FOR_EACH(entry, IntDynamicVector) if(entry == 0) entry++;
  FOR_EACH(entry, RealStaticVector) if(isEqual(entry, Zero)) entry++;
  FOR_EACH(entry, RealDynamicVector) if(isEqual(entry, Zero)) entry++;

  IntStaticVectorTest = IntStaticVector / IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector / IntDynamicVector;
   RealStaticVectorTest = RealStaticVector / RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector / RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 1);
    EXPECT_EQ(IntDynamicVectorTest[i], 1);
    EXPECT_EQ(RealStaticVectorTest[i], One);
    EXPECT_EQ(RealDynamicVectorTest[i], One);
  }
}

TEST_F(VectorTest, VectorPlusEqual)
{
  IntStaticVectorTest = IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector;
   RealStaticVectorTest = RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector;

  IntStaticVectorTest += IntStaticVector;
  IntDynamicVectorTest += IntDynamicVector;
   RealStaticVectorTest += RealStaticVector;
   RealDynamicVectorTest += RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 2 * IntStaticVector[i]);
    EXPECT_EQ(IntDynamicVectorTest[i], 2 * IntDynamicVector[i]);
    EXPECT_EQ(RealStaticVectorTest[i], Two * RealStaticVector[i]);
    EXPECT_EQ(RealDynamicVectorTest[i], Two * RealDynamicVector[i]);
  }
}

TEST_F(VectorTest, VectorMinusEqual)
{
  IntStaticVectorTest = IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector;
   RealStaticVectorTest = RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector;

  IntStaticVectorTest -= IntStaticVector;
  IntDynamicVectorTest -= IntDynamicVector;
   RealStaticVectorTest -= RealStaticVector;
   RealDynamicVectorTest -= RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 0);
    EXPECT_EQ(IntDynamicVectorTest[i], 0);
    EXPECT_EQ(RealStaticVectorTest[i], Zero);
    EXPECT_EQ(RealDynamicVectorTest[i], Zero);
  }
}

TEST_F(VectorTest, VectorMultiplyEqual)
{
  IntStaticVectorTest = IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector;
   RealStaticVectorTest = RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector;

  IntStaticVectorTest *= IntStaticVector;
  IntDynamicVectorTest *= IntDynamicVector;
   RealStaticVectorTest *= RealStaticVector;
   RealDynamicVectorTest *= RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], iPow(IntStaticVector[i], 2));
    EXPECT_EQ(IntDynamicVectorTest[i], iPow(IntDynamicVector[i], 2));
    EXPECT_EQ(RealStaticVectorTest[i], iPow(RealStaticVector[i], 2));
    EXPECT_EQ(RealDynamicVectorTest[i], iPow(RealDynamicVector[i], 2));
  }
}

TEST_F(VectorTest, VectorDivideEqual)
{
  FOR_EACH(entry, IntStaticVector) if(entry == 0) entry++;
  FOR_EACH(entry, IntDynamicVector) if(entry == 0) entry++;
  FOR_EACH(entry, RealStaticVector) if(isEqual(entry, Zero)) entry++;
  FOR_EACH(entry, RealDynamicVector) if(isEqual(entry, Zero)) entry++;

  IntStaticVectorTest = IntStaticVector;
  IntDynamicVectorTest = IntDynamicVector;
   RealStaticVectorTest = RealStaticVector;
   RealDynamicVectorTest = RealDynamicVector;

  IntStaticVectorTest /= IntStaticVector;
  IntDynamicVectorTest /= IntDynamicVector;
   RealStaticVectorTest /= RealStaticVector;
   RealDynamicVectorTest /= RealDynamicVector;

  FOR(i, ContainerSize)
  {
    EXPECT_EQ(IntStaticVectorTest[i], 1);
    EXPECT_EQ(IntDynamicVectorTest[i], 1);
    EXPECT_EQ(RealStaticVectorTest[i], One);
    EXPECT_EQ(RealDynamicVectorTest[i], One);
  }
}

/***************************************************************************************************************************************************************
* Test Other Vector Operations
***************************************************************************************************************************************************************/
TEST_F(VectorTest, InnerProduct)
{
  // Parallel test (same/opposite directions)
  RealStaticVector = Normalise(RealStaticVector);
  EXPECT_NEAR(InnerProduct(RealStaticVector, RealStaticVector), One, Small);
  EXPECT_NEAR(InnerProduct(RealStaticVector, -RealStaticVector), -One, Small);

  // Orthogonal test
  EXPECT_DOUBLE_EQ(InnerProduct(xAxis2, yAxis2), Zero);
  EXPECT_DOUBLE_EQ(InnerProduct(xAxis3, yAxis3), Zero);
  EXPECT_DOUBLE_EQ(InnerProduct(xAxis3, zAxis3), Zero);
  EXPECT_DOUBLE_EQ(InnerProduct(yAxis3, zAxis3), Zero);

  // Arbitrary test
  SVectorR3 random0, random1;
  random0.Randomise();
  random1.Randomise();
  EXPECT_NEAR(InnerProduct(random0, random1), L2Norm(random0) * L2Norm(random1) * std::cos(ComputeAngle(random0, random1)), Two*Small);
}

TEST_F(VectorTest, CrossProduct)
{
  SVectorR3 random0;
  random0.Randomise();
  SVectorR3 cross_product;

  // Parallel test (same/opposite directions)
  cross_product = CrossProduct(ToVector<2>(random0), ToVector<2>(random0));
  FOR_EACH(entry, cross_product) EXPECT_DOUBLE_EQ(entry, Zero);
  cross_product = CrossProduct(ToVector<2>(random0), ToVector<2>(-random0));
  FOR_EACH(entry, cross_product) EXPECT_DOUBLE_EQ(entry, Zero);
  cross_product = CrossProduct(random0, random0);
  FOR_EACH(entry, cross_product) EXPECT_DOUBLE_EQ(entry, Zero);
  cross_product = CrossProduct(random0, -random0);
  FOR_EACH(entry, cross_product) EXPECT_DOUBLE_EQ(entry, Zero);

  // Orthogonal test
  cross_product = CrossProduct(xAxis2, yAxis2);
  EXPECT_TRUE(cross_product == zAxis3);
  cross_product = CrossProduct(xAxis3, yAxis3);
  EXPECT_TRUE(cross_product == zAxis3);
  cross_product = CrossProduct(yAxis3, xAxis3);
  EXPECT_TRUE(cross_product == -zAxis3);
  cross_product = CrossProduct(yAxis3, zAxis3);
  EXPECT_TRUE(cross_product == xAxis3);
  cross_product = CrossProduct(zAxis3, xAxis3);
  EXPECT_TRUE(cross_product == yAxis3);

  // Arbitrary test
  SVectorR3 random1;
  random1.Randomise();
  EXPECT_NEAR(Magnitude(CrossProduct(random0, random1)), Magnitude(random0) * Magnitude(random1) * std::sin(ComputeAngle(random0, random1)), Two*Small);
}

TEST_F(VectorTest, L1Norm)
{
  EXPECT_DOUBLE_EQ(L1Norm(RealStaticVector), std::accumulate(RealStaticVector.begin(), RealStaticVector.end(), Zero));
  EXPECT_EQ(L1Norm(IntStaticVector), std::accumulate(IntStaticVector.begin(), IntStaticVector.end(), 0));
}

TEST_F(VectorTest, L2Norm)
{
  Real inner_prod{};
  FOR_EACH(entry, RealStaticVector) inner_prod += iPow(entry, 2);
  EXPECT_DOUBLE_EQ(L2Norm(RealStaticVector), std::sqrt(inner_prod));

  inner_prod = Zero;
  FOR_EACH(entry, IntStaticVector) inner_prod += iPow(entry, 2);
  EXPECT_DOUBLE_EQ(L2Norm(IntStaticVector), std::sqrt(inner_prod));
}

TEST_F(VectorTest, LInfNorm)
{
  EXPECT_EQ(LInfNorm(RealStaticVector), *std::max_element(RealStaticVector.begin(), RealStaticVector.end()));
  EXPECT_EQ(LInfNorm(IntStaticVector), *std::max_element(IntStaticVector.begin(), IntStaticVector.end()));
}

TEST_F(VectorTest, isNormalised)
{
  EXPECT_FALSE(isNormalised(RealStaticVector));
   RealStaticVector = Normalise(RealStaticVector);
  EXPECT_TRUE(isNormalised(RealStaticVector));
}

TEST_F(VectorTest, Normalise)
{
   RealStaticVector = Normalise(RealStaticVector);
  Real inner_prod{};
  FOR_EACH(entry, RealStaticVector) inner_prod += iPow(entry, 2);
  EXPECT_DOUBLE_EQ(std::sqrt(inner_prod), One);
}

TEST_F(VectorTest, ComputeAngle)
{
  // Non-oriented tests
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis2, xAxis2), Zero);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis2, -xAxis2), Pi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis2, SVectorR2{One, One}), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(SVectorR2{One, One}, xAxis2), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis2, SVectorR2{-One, One}), ThreeQuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(SVectorR2{-One, One}, xAxis2), ThreeQuarterPi);

  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, xAxis3), Zero);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, -xAxis3), Pi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, SVectorR3{One, One, Zero}), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, SVectorR3{-One, One, Zero}), ThreeQuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, yAxis3), HalfPi);
  EXPECT_DOUBLE_EQ(ComputeAngle(xAxis3, -yAxis3), HalfPi);

  // Oriented tests
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, xAxis2), Zero);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, -xAxis2), Pi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, SVectorR2{One, One}), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, SVectorR2{-One, One}), ThreeQuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, SVectorR2{-One, One}, -zAxis3), -ThreeQuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis2, SVectorR2{-One, -One}), -ThreeQuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(SVectorR2{One, One}, xAxis2), -QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(SVectorR2{One, One}, xAxis2, -zAxis3), QuarterPi);

  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis3, xAxis3), Zero);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis3, -xAxis3), Pi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis3, SVectorR3{One, One, Zero}), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(SVectorR3{One, One, Zero}, xAxis3), -QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(SVectorR3{One, One, Zero}, xAxis3, -zAxis3), QuarterPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis3, yAxis3), HalfPi);
  EXPECT_DOUBLE_EQ(ComputeAngle<true>(xAxis3, -yAxis3), -HalfPi);
}

TEST_F(VectorTest, isAligned)
{
  EXPECT_TRUE(isAligned(xAxis2, xAxis2));
  EXPECT_TRUE(isAligned(xAxis2, -xAxis2));
  EXPECT_FALSE(isAligned(xAxis2, yAxis2));
  EXPECT_FALSE(isAligned(xAxis2, -yAxis2));
  EXPECT_FALSE(isAligned(xAxis2, SVectorR2{One, One}));
  EXPECT_TRUE(isAligned(xAxis2, SVectorR2{One, One}, DegToRad(45.00001)));
  EXPECT_TRUE(isAligned(xAxis2, SVectorR2{Ten, One}));

  EXPECT_TRUE(isAligned(xAxis3, xAxis3));
  EXPECT_TRUE(isAligned(xAxis3, -xAxis3));
  EXPECT_FALSE(isAligned(xAxis3, yAxis3));
  EXPECT_FALSE(isAligned(xAxis3, -yAxis3));
  EXPECT_FALSE(isAligned(xAxis3, SVectorR3{One, One, Zero}));
  EXPECT_TRUE(isAligned(xAxis3, SVectorR3{One, One, Zero}, DegToRad(45.00001)));
  EXPECT_TRUE(isAligned(xAxis3, SVectorR3{Ten, One, Zero}));

  EXPECT_THROW(isAligned(xAxis3, SVectorR3{One, Zero, Zero}, DegToRad(90.00001)), std::domain_error);
}

}

#endif