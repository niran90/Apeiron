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

#include <gtest/gtest.h>
#include "../include/Global.h"
#include "../include/Random.h"

namespace aprn{

class ApeironTest : public testing::Test
{
public:
  Random<int>  RandomInt;
  Random<Real> RandomReal;

  std::array<int, 100>  IntArray;
  std::array<Real, 100> RealArray;

  ApeironTest()
    : RandomInt(-10, 10), RandomReal(-Ten, Ten) {}

  void SetUp() override
  {
    FOR_EACH(entry, IntArray)  entry = RandomInt();
    FOR_EACH(entry, RealArray) entry = RandomReal();
  }
};

}
