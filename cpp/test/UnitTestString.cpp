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

TEST_F(ApeironTest, ToString)
{
   std::string str = ToString(1);
   EXPECT_EQ(str, "1");

   str = ToString(5) + "s";
   EXPECT_EQ(str, "5s");

   str = ToString(10.5) + "m";
   EXPECT_EQ(str, "10.5m");
}

TEST_F(ApeironTest, ToNumber)
{
   EXPECT_EQ(ToNumber("1")           , 1);
   EXPECT_EQ(ToNumber<int>("1")      , 1);
   EXPECT_EQ(ToNumber<long>("1")     , 1);
   EXPECT_EQ(ToNumber<long long>("1"), 1);

   EXPECT_EQ(ToNumber<float>("5.5")      , 5.5);
   EXPECT_EQ(ToNumber<double>("5.5")     , 5.5);
   EXPECT_EQ(ToNumber<long double>("5.5"), 5.5);
   EXPECT_EQ(ToNumber<double>("5.5e-3")  , 0.0055);
   EXPECT_EQ(ToNumber<double>("5.5e10")  , 5.5e10);
}

TEST_F(ApeironTest, isSubstring)
{
   EXPECT_TRUE(isSubstring("T"    , "True"));
   EXPECT_TRUE(isSubstring("r"    , "True"));
   EXPECT_TRUE(isSubstring("u"    , "True"));
   EXPECT_TRUE(isSubstring("e"    , "True"));
   EXPECT_TRUE(isSubstring("rue"  , "True"));
   EXPECT_TRUE(isSubstring(""     , "True"));
   EXPECT_TRUE(isSubstring("wo wo", "Two words"));
   EXPECT_TRUE(isSubstring(" "    , "Two words"));
   EXPECT_TRUE(isSubstring("\t"   , "Two\twords"));
   EXPECT_TRUE(isSubstring("\n"   , "Two\nwords"));

   EXPECT_FALSE(isSubstring("T"   , "False"));
   EXPECT_FALSE(isSubstring("r"   , "False"));
   EXPECT_FALSE(isSubstring("u"   , "False"));
   EXPECT_FALSE(isSubstring("m"   , "False"));
   EXPECT_FALSE(isSubstring("p"   , "False"));
   EXPECT_FALSE(isSubstring("Pal" , "False"));
   EXPECT_FALSE(isSubstring("wo o", "Two words"));
   EXPECT_FALSE(isSubstring("\t"  , "Two words"));
   EXPECT_FALSE(isSubstring("\n"  , "Two words"));
   EXPECT_FALSE(isSubstring(" "   , "Two\twords"));
   EXPECT_FALSE(isSubstring("\n"  , "Two\twords"));
   EXPECT_FALSE(isSubstring(" "   , "Two\nwords"));
   EXPECT_FALSE(isSubstring("\t"  , "Two\nwords"));

   // Test raw string literals
   EXPECT_TRUE(isSubstring("\\pi"  , R"(2\pi)"));
   EXPECT_TRUE(isSubstring(R"(\pi)", R"(2\pi)"));
   EXPECT_TRUE(isSubstring(R"(\pi)", "2\\pi"));
   EXPECT_TRUE(isSubstring(R"(\)"  , R"(2\pi)"));
   EXPECT_TRUE(isSubstring("\\"    , R"(2\pi)"));
   EXPECT_TRUE(isSubstring("   "   , R"(Two   words)"));
   EXPECT_TRUE(isSubstring("\n    ", R"(Two
    words)"));
}

TEST_F(ApeironTest, Replace)
{
   EXPECT_EQ(Replace("T"    , "D" , "True")     , "Drue");
   EXPECT_EQ(Replace("r"    , "s" , "True")     , "Tsue");
   EXPECT_EQ(Replace("u"    , "e" , "True")     , "Tree");
   EXPECT_EQ(Replace("e"    , "D" , "True")     , "TruD");
   EXPECT_EQ(Replace("rue"  , "oe", "True")     , "Toe");
   EXPECT_EQ(Replace(""     , "-" , "True")     , "True");
   EXPECT_EQ(Replace("wo wo", "u" , "Two words"), "Turds");
}

TEST_F(ApeironTest, Remove)
{
   EXPECT_EQ(Remove("T"    , "True")     , "rue");
   EXPECT_EQ(Remove("r"    , "True")     , "Tue");
   EXPECT_EQ(Remove("u"    , "True")     , "Tre");
   EXPECT_EQ(Remove("e"    , "True")     , "Tru");
   EXPECT_EQ(Remove("rue"  , "True")     , "T");
   EXPECT_EQ(Remove(""     , "True")     , "True");
   EXPECT_EQ(Remove("wo wo", "Two words"), "Trds");
   EXPECT_EQ(Remove("{"    , "args{1}{2}"), "args1}2}");
}

TEST_F(ApeironTest, Split)
{
   std::vector<std::string> sub_strings;
   sub_strings.push_back("This");
   sub_strings.push_back("is");
   sub_strings.push_back("the");
   sub_strings.push_back("string");
   sub_strings.push_back("that");
   sub_strings.push_back("I");
   sub_strings.push_back("would");
   sub_strings.push_back("like");
   sub_strings.push_back("to");
   sub_strings.push_back("test.");

   std::vector<std::string> delimiters;
   delimiters.push_back(",");
   delimiters.push_back("-");
   delimiters.push_back("_");
   delimiters.push_back("(");
   delimiters.push_back("{");
   delimiters.push_back("[");

   const auto n_words = sub_strings.size();
   FOR_EACH(delimiter, delimiters)
   {
      std::string test_string = sub_strings.front();
      FOR(i, 1, n_words) test_string = test_string + delimiter + sub_strings[i];

      std::vector<std::string> test_sub_strings = Split(test_string, delimiter);
      EXPECT_EQ(test_sub_strings.size(), n_words);
      FOR(i, n_words) EXPECT_EQ(test_sub_strings[i], sub_strings[i]);
   }

   EXPECT_EQ(Split("Test character split.").size(), 21);
}

TEST_F(ApeironTest, GetFirstEnclosure)
{
   std::string enclosure;
   bool was_found;

   std::tie(enclosure, was_found) = GetFirstEnclosure("Test (this).", '(', ')');
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "this");
   std::tie(enclosure, was_found) = GetFirstEnclosure("Test (this).", '(', ')', true);
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "(this)");

   std::tie(enclosure, was_found) = GetFirstEnclosure("Now [test] (this).", '[', ']');
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "test");
   std::tie(enclosure, was_found) = GetFirstEnclosure("Now [test] (this).", '[', ']', true);
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "[test]");

   std::tie(enclosure, was_found) = GetFirstEnclosure("Now (test (this) nested) case.", '(', ')');
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "test (this) nested");
   std::tie(enclosure, was_found) = GetFirstEnclosure("Now (test (this) nested) case.", '(', ')', true);
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "(test (this) nested)");

   std::tie(enclosure, was_found) = GetFirstEnclosure("Test $this$ too.", '$', '$');
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "this");
   std::tie(enclosure, was_found) = GetFirstEnclosure("Test $this$ too.", '$', '$', true);
   EXPECT_TRUE(was_found);
   EXPECT_EQ(enclosure, "$this$");

   std::tie(enclosure, was_found) = GetFirstEnclosure("Next, test (this).", '[', ']');
   EXPECT_FALSE(was_found);
   EXPECT_EQ(enclosure, "");
   std::tie(enclosure, was_found) = GetFirstEnclosure("Next, test (this).", '[', ']', true);
   EXPECT_FALSE(was_found);
   EXPECT_EQ(enclosure, "");


   EXPECT_DEATH(GetFirstEnclosure("Also test $this.", '$', '$'), "");
   EXPECT_DEATH(GetFirstEnclosure("Also test $this.", '$', '$', true), "");
   EXPECT_DEATH(GetFirstEnclosure("Finally, test {this.", '{', '}'), "");
   EXPECT_DEATH(GetFirstEnclosure("Finally, test {this.", '{', '}', true), "");
}

TEST_F(ApeironTest, GetFirstEnclosureChain)
{
   std::vector<std::string> enclosure_chain;

   enclosure_chain = GetFirstEnclosureChain("Test (this) (first).", '(', ')');
   EXPECT_EQ(enclosure_chain.size(), 1);
   EXPECT_EQ(enclosure_chain[0], "this");
   enclosure_chain = GetFirstEnclosureChain("Test (this) (first).", '(', ')', true);
   EXPECT_EQ(enclosure_chain.size(), 1);
   EXPECT_EQ(enclosure_chain[0], "(this)");

   enclosure_chain = GetFirstEnclosureChain("Test [this][sentence] [next].", '[', ']');
   EXPECT_EQ(enclosure_chain.size(), 2);
   EXPECT_EQ(enclosure_chain[0], "this");
   EXPECT_EQ(enclosure_chain[1], "sentence");
   enclosure_chain = GetFirstEnclosureChain("Test [this][sentence] [next].", '[', ']', true);
   EXPECT_EQ(enclosure_chain.size(), 2);
   EXPECT_EQ(enclosure_chain[0], "[this]");
   EXPECT_EQ(enclosure_chain[1], "[sentence]");

   enclosure_chain = GetFirstEnclosureChain("[test [this] [nested]][case] [also].", '[', ']');
   EXPECT_EQ(enclosure_chain.size(), 2);
   EXPECT_EQ(enclosure_chain[0], "test [this] [nested]");
   EXPECT_EQ(enclosure_chain[1], "case");
   enclosure_chain = GetFirstEnclosureChain("[test [this] [nested]][case] [also].", '[', ']', true);
   EXPECT_EQ(enclosure_chain.size(), 2);
   EXPECT_EQ(enclosure_chain[0], "[test [this] [nested]]");
   EXPECT_EQ(enclosure_chain[1], "[case]");

   enclosure_chain = GetAllEnclosures("Next, (test)(this).", '[', ']');
   EXPECT_TRUE(enclosure_chain.empty());
   enclosure_chain = GetAllEnclosures("Next, (test)(this).", '[', ']', true);
   EXPECT_TRUE(enclosure_chain.empty());

   EXPECT_DEATH(GetAllEnclosures("Also {test}{this.", '{', '}'), "");
   EXPECT_DEATH(GetAllEnclosures("Also {test}{this.", '{', '}', true), "");
}

TEST_F(ApeironTest, GetAllEnclosures)
{
   std::vector<std::string> enclosures;

   enclosures = GetAllEnclosures("Test (this) (first).", '(', ')');
   EXPECT_EQ(enclosures.size(), 2);
   EXPECT_EQ(enclosures[0], "this");
   EXPECT_EQ(enclosures[1], "first");
   enclosures = GetAllEnclosures("Test (this) (first).", '(', ')', true);
   EXPECT_EQ(enclosures.size(), 2);
   EXPECT_EQ(enclosures[0], "(this)");
   EXPECT_EQ(enclosures[1], "(first)");

   enclosures = GetAllEnclosures("[Test] [this] [second].", '[', ']');
   EXPECT_EQ(enclosures.size(), 3);
   EXPECT_EQ(enclosures[0], "Test");
   EXPECT_EQ(enclosures[1], "this");
   EXPECT_EQ(enclosures[2], "second");
   enclosures = GetAllEnclosures("[Test] [this] [second].", '[', ']', true);
   EXPECT_EQ(enclosures.size(), 3);
   EXPECT_EQ(enclosures[0], "[Test]");
   EXPECT_EQ(enclosures[1], "[this]");
   EXPECT_EQ(enclosures[2], "[second]");

   enclosures = GetAllEnclosures("[Now], [test [this] [nested]] [case].", '[', ']');
   EXPECT_EQ(enclosures.size(), 3);
   EXPECT_EQ(enclosures[0], "Now");
   EXPECT_EQ(enclosures[1], "test [this] [nested]");
   EXPECT_EQ(enclosures[2], "case");
   enclosures = GetAllEnclosures("[Now], [test [this] [nested]] [case].", '[', ']', true);
   EXPECT_EQ(enclosures.size(), 3);
   EXPECT_EQ(enclosures[0], "[Now]");
   EXPECT_EQ(enclosures[1], "[test [this] [nested]]");
   EXPECT_EQ(enclosures[2], "[case]");

   enclosures = GetAllEnclosures("Next, test (this).", '[', ']');
   EXPECT_TRUE(enclosures.empty());
   enclosures = GetAllEnclosures("Next, test (this).", '[', ']', true);
   EXPECT_TRUE(enclosures.empty());

   EXPECT_DEATH(GetAllEnclosures("Also test $this.", '$', '$'), "");
   EXPECT_DEATH(GetAllEnclosures("Also test $this.", '$', '$', true), "");
   EXPECT_DEATH(GetAllEnclosures("Finally, test {this.", '{', '}'), "");
   EXPECT_DEATH(GetAllEnclosures("Finally, test {this.", '{', '}', true), "");
}

}

#endif