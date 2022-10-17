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
#include "../include/Curve.h"

#ifdef DEBUG_MODE

namespace aprn::mnfld {

/***************************************************************************************************************************************************************
* Curve Test Fixture
***************************************************************************************************************************************************************/
class CurveTest : public testing::Test
{
public:
  Random<int>  RandomInt;
  Random<Real> RandomReal;

  CurveTest()
    : RandomInt(-10, 10), RandomReal(-Ten, Ten) {}

  void
  SetUp() override
  {

  }
};


/***************************************************************************************************************************************************************
* Linear/Piecewise Linear Curves
***************************************************************************************************************************************************************/
TEST_F(CurveTest, Line)
{
  SVectorR3 p, direction, norm_direction, centre;
  direction.Randomise();
  centre.Randomise();
  norm_direction = Normalise(direction);
  Line line(direction, centre);

  // Not unit speed parametrised
  p = line.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], centre[i]);

  p = line.Point(One);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], centre[i] + direction[i]);

  p = line.Point(-One);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], centre[i] - direction[i]);

  // Unit speed parametrised
  line.SetIfUnitSpeed(true);
  const Real random = RandomReal();
  p = line.Point(random);
  FOR(i, 3) EXPECT_NEAR(p[i], centre[i] + random * norm_direction[i], Two * Small);

  p = line.Point(-random);
  FOR(i, 3) EXPECT_NEAR(p[i], centre[i] - random * norm_direction[i], Two * Small);
}

TEST_F(CurveTest, Ray)
{
  SVectorR3 p, start, direction, norm_direction;
  start.Randomise();
  direction.Randomise();
  norm_direction = Normalise(direction);
  Ray ray(direction, start);

  // Not unit speed parametrised
  p = ray.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], start[i]);

  p = ray.Point(One);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], start[i] + direction[i]);

  p = ray.Point(Two);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], start[i] + Two * direction[i]);

  EXPECT_THROW(ray.Point(-Small), std::domain_error);

  // Unit speed parametrised
  ray.SetIfUnitSpeed(true);
  RandomReal.Reset(Zero, Ten);
  const Real random = RandomReal();
  p = ray.Point(random);
  FOR(i, 3) EXPECT_NEAR(p[i], start[i] + random * norm_direction[i], Small);

  EXPECT_THROW(ray.Point(-Small), std::domain_error);
}

TEST_F(CurveTest, Segment)
{
  SVectorR3 p, start, end, direction, norm_direction;
  start.Randomise();
  end.Randomise();
  direction = end - start;
  norm_direction = Normalise(direction);
  const Real magnitude = Magnitude(direction);
  Segment segment(start, end);

  // Not unit speed parametrised
  p = segment.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], start[i]);

  p = segment.Point(One);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], end[i]);

  p = segment.Point(Half);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], start[i] + Half * direction[i]);

  EXPECT_THROW(segment.Point(-Ten * Small), std::domain_error);
  EXPECT_THROW(segment.Point(One + Ten * Small), std::domain_error);

  // Unit speed parametrised
  segment.SetIfUnitSpeed(true);
  RandomReal.Reset(Zero, magnitude);
  const Real random = RandomReal();
  p = segment.Point(random);
  FOR(i, 3) EXPECT_NEAR(p[i], start[i] + random * norm_direction[i], Small);

  EXPECT_THROW(segment.Point(-Ten * Small), std::domain_error);
  EXPECT_THROW(segment.Point(magnitude + Ten * Small), std::domain_error);
}

TEST_F(CurveTest, SegmentChain)
{
  RandomInt.Reset(3, 12);
  const size_t n_vertices = RandomInt();

  DynamicArray<SVectorR3> vertices(n_vertices);
  DynamicArray<SVectorR3> directions(n_vertices - 1);
  DynamicArray<SVectorR3> norm_directions(n_vertices - 1);
  DynamicArray<Real> magnitudes(n_vertices - 1);
  Real chain_length{};
  FOR(i, vertices.size())
  {
    vertices[i].Randomise();
    vertices[i] *= Ten;
    if(i > 0)
    {
      directions[i - 1] = vertices[i] - vertices[i - 1];
      norm_directions[i - 1] = Normalise(directions[i - 1]);
      magnitudes[i - 1] = Magnitude(directions[i - 1]);
      chain_length += magnitudes[i - 1];
    }
  }

  SegmentChain chain(vertices);
  SVectorR3 p;

  // Not unit speed parametrised
  p = chain.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], vertices.front()[i]);

  p = chain.Point(One);
  FOR(i, 3) EXPECT_NEAR(p[i], vertices.back()[i], Ten * Small);

  RandomReal.Reset(Zero, One);
  Real random = RandomReal();
  p = chain.Point(random);
  Real param_length = random * chain_length;
  Real sum{};
  SVectorR3 p_check;
  bool isfound{};
  FOR(i, magnitudes.size())
  {
    const auto sum_prev = sum;
    sum += magnitudes[i];
    if(param_length <= sum)
    {
      const Real r = (param_length - sum_prev) / magnitudes[i];
      p_check = (One - r) * vertices[i] + r * vertices[i + 1];
      isfound = true;
      break;
    }
  }
  EXPECT_TRUE(isfound);
  FOR(i, 3) EXPECT_NEAR(p[i], p_check[i], Two * Small);

  EXPECT_THROW(chain.Point(-Small), std::domain_error);
  EXPECT_THROW(chain.Point(One + Small), std::domain_error);

  // Unit speed parametrised
  chain.SetIfUnitSpeed(true);
  RandomReal.Reset(Zero, chain_length);
  random = RandomReal();
  p = chain.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], vertices.front()[i]);

  p = chain.Point(chain_length);
  FOR(i, 3) EXPECT_NEAR(p[i], vertices.back()[i], Ten * Small);

  RandomReal.Reset(Zero, chain_length);
  random = RandomReal();
  p = chain.Point(random);
  param_length = random;
  sum = Zero;
  isfound = false;
  FOR(i, magnitudes.size())
  {
    const auto sum_prev = sum;
    sum += magnitudes[i];
    if(param_length <= sum)
    {
      const Real r = (param_length - sum_prev) / magnitudes[i];
      p_check = (One - r) * vertices[i] + r * vertices[i + 1];
      isfound = true;
      break;
    }
  }
  EXPECT_TRUE(isfound);
  FOR(i, 3) EXPECT_NEAR(p[i], p_check[i], Two * Small);

  EXPECT_THROW(chain.Point(-Small), std::domain_error);
  EXPECT_THROW(chain.Point(chain_length + Ten * Small), std::domain_error);

  // Test closed chain
  chain = SegmentChain(vertices, true);
  p = chain.Point(Zero);
  FOR(i, 3) EXPECT_DOUBLE_EQ(p[i], vertices.front()[i]);

  p = chain.Point(One);
  FOR(i, 3) EXPECT_NEAR(p[i], vertices.front()[i], Ten * Small);
}

/***************************************************************************************************************************************************************
* Circular/Elliptical Curves
***************************************************************************************************************************************************************/
TEST_F(CurveTest, Circle)
{
  RandomReal.Reset(One, Ten);
  const Real radius = RandomReal();
  SVectorR2 p, p_check, centre;
  centre.Randomise();
  Circle circle(radius, centre);

  // Not unit speed parametrised
  p = circle.Point(Zero);
  FOR(i, 2) EXPECT_DOUBLE_EQ(p[i], centre[i] + radius * xAxis2[i]);

  p = circle.Point(QuartPi);
  FOR(i, 2) EXPECT_NEAR(p[i], centre[i] + radius / Sqrt(Two), Two * Small);

  p = circle.Point(HalfPi);
  FOR(i, 2) EXPECT_NEAR(p[i], centre[i] + radius * yAxis2[i], Two * Small);

  p = circle.Point(-HalfPi);
  FOR(i, 2) EXPECT_DOUBLE_EQ(p[i], centre[i] - radius * yAxis2[i]);

  // Unit speed parametrised
  circle.SetIfUnitSpeed(true);
  RandomReal.Reset(-Ten, Ten);
  const Real random = RandomReal();
  p = circle.Point(random);
  Real theta = random / radius;
  p_check = centre + radius * SVectorR2{Cos(theta), Sin(theta)};
  FOR(i, 2) EXPECT_NEAR(p[i], p_check[i], 30.0 * Small);

  p = circle.Point(-random);
  p_check = centre + radius * SVectorR2{Cos(theta), -Sin(theta)};
  FOR(i, 2) EXPECT_NEAR(p[i], p_check[i], 30.0 * Small);
}

TEST_F(CurveTest, Ellipse)
{
  RandomReal.Reset(One, Ten);
  const Real radius_x = RandomReal();
  const Real radius_y = RandomReal();
  SVectorR2 p, p_check, centre;
  centre.Randomise();
  Ellipse ellipse(radius_x, radius_y, centre);

  // Not unit speed parametrised
  p = ellipse.Point(Zero);
  FOR(i, 2) EXPECT_DOUBLE_EQ(p[i], centre[i] + radius_x * xAxis2[i]);

  p = ellipse.Point(QuartPi);
  p_check = centre + SVectorR2{radius_x / Sqrt(Two), radius_y / Sqrt(Two)};
  FOR(i, 2) EXPECT_NEAR(p[i], p_check[i], Two * Small);

  p = ellipse.Point(HalfPi);
  FOR(i, 2) EXPECT_NEAR(p[i], centre[i] + radius_y * yAxis2[i], Two * Small);

  p = ellipse.Point(-HalfPi);
  FOR(i, 2) EXPECT_NEAR(p[i], centre[i] - radius_y * yAxis2[i], Two * Small);

  // Unit speed parametrised - requires root-finding and quadrature first.
}

}

#endif