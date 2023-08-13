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

#include "../../../include/Global.h"
#include "DataContainer/include/Array.h"

#include <filesystem>
#include <stdlib.h>
#include <string_view>

namespace aprn::flmgr {

namespace fs = std::filesystem;
typedef fs::path Path;

/***************************************************************************************************************************************************************
* File Handling Functions
***************************************************************************************************************************************************************/
Path FileName(const std::string_view& file_path);

bool FileExists(const Path& file_path);

bool FileIsEmpty(const Path& file_path);

void ClearFile(const Path& file_path);

bool DeleteFile(const Path& file_path);

void MoveFile(const Path& from_file_path, const Path& to_dir_path);

void CopyFile(const Path& from_file_path, const Path& to_dir_path);

/***************************************************************************************************************************************************************
* Directory Handling Functions
***************************************************************************************************************************************************************/
bool isDirectory(const Path& path);

bool DirectoryExists(const Path& dir_path);

void ClearDirectory(const Path& dir_path);

void CreateDirectory(const Path& dir_path, bool clear_if_exists = false);

void CreateDirectories(const Path& dir_path);

bool DirectoryIsEmpty(const Path& dir_path);

bool DeleteDirectory(const Path& dir_path);

void CopyDirectory(const Path& from_dir_path, const Path& to_dir_path);

void MoveDirectory(const std::string_view& from_dir_path, const std::string_view& to_dir_path);

/***************************************************************************************************************************************************************
* General Linux Commands
***************************************************************************************************************************************************************/
int RunCommand(const std::string& cmd, bool show_output = false);

int RunCommandFrom(const std::string& cmd, const Path& dir_path, bool show_output = false);

void RunCommands(const DArray<std::string>& cmds, bool show_output = false);

void RunCommandsFrom(const DArray<std::string>& cmds, const Path& dir_path, bool show_output = false);

bool CommandExists(const std::string& cmd);

void CompileTeXFile(const std::string& tex_compiler, const Path& tex_path, bool show_output = false);

void ConvertPDFtoPNG(Path pdf_path, const unsigned pixel_density = 600, bool show_output = false);

}
