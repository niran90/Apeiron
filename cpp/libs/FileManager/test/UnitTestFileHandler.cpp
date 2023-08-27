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
#include "../include/FileSystem.h"

#include <filesystem>
#include <string_view>

#ifdef DEBUG_MODE

namespace aprn::flmgr {

namespace fs = std::filesystem;

class FileHandlerTest : public testing::Test
{
 public:
   const std::string DataDir{"../libs/FileManager/test/data"}; // Relative to CMake build directory.

   FileHandlerTest() {}

   void
   SetUp() override {}
};

/***************************************************************************************************************************************************************
* Test File Handling Functions
***************************************************************************************************************************************************************/
TEST_F(FileHandlerTest, FileExists)
{
   EXPECT_TRUE(FileExists(DataDir + "/empty_file.txt"));
   EXPECT_FALSE(FileExists(DataDir + "/image.png"));
}

TEST_F(FileHandlerTest, FileIsEmpty)
{
   EXPECT_TRUE(FileIsEmpty(DataDir + "/empty_file.txt"));
   EXPECT_FALSE(FileIsEmpty(DataDir + "/non_empty_directory/non_empty_file.txt"));
   EXPECT_DEATH(FileIsEmpty(DataDir + "/image.png"), "");
}

TEST_F(FileHandlerTest, CopyFile)
{

}

TEST_F(FileHandlerTest, MoveFile)
{

}

/***************************************************************************************************************************************************************
* Test Directory Handling Functions
***************************************************************************************************************************************************************/
TEST_F(FileHandlerTest, isDirectory)
{
   EXPECT_TRUE(isDirectory(DataDir + "/empty_directory")); // Existing directory
   EXPECT_FALSE(isDirectory(DataDir + "/empty_file.txt")); // Existing file
   EXPECT_FALSE(isDirectory(DataDir + "/images"));         // Non-existant file/directory
}

TEST_F(FileHandlerTest, DirectoryExists)
{
   EXPECT_TRUE(DirectoryExists(DataDir + "/empty_directory"));
   EXPECT_FALSE(DirectoryExists(DataDir + "/data"));
}

TEST_F(FileHandlerTest, DirectoryIsEmpty)
{
   EXPECT_FALSE(DirectoryIsEmpty(DataDir + "/non_empty_directory"));
   EXPECT_DEATH(DirectoryIsEmpty(DataDir + "/non_existant_directory"), "");
}

TEST_F(FileHandlerTest, CopyDirectory)
{

}

TEST_F(FileHandlerTest, MoveDirectory)
{

}

}

#endif