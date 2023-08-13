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

#include "../include/FileSystem.h"

namespace aprn::flmgr {

/***************************************************************************************************************************************************************
* File Handling Functions
***************************************************************************************************************************************************************/
Path
FileName(const std::string_view& file_path) { return Path(file_path).filename(); }

bool
FileExists(const Path& file_path) { return fs::exists(file_path); }

bool
FileIsEmpty(const Path& file_path)
{
   ASSERT(FileExists(file_path), "The file/directory ", file_path.filename(), " does not exist.")
   return fs::is_empty(file_path);
}

void
ClearFile(const Path& file_path)
{
   ASSERT(FileExists(file_path), "The file ", file_path.filename(), " does not exist.")
   fs::resize_file(file_path, 0);
}

bool
DeleteFile(const Path& file_path)
{
   ASSERT(FileExists(file_path), "The file/directory ", file_path.filename(), " does not exist.")
   return fs::remove(file_path);
}

void
MoveFile(const Path& from_file_path, const Path& to_dir_path)
{
   ASSERT(FileExists(from_file_path), "The file ", from_file_path.filename(), " does not exist.")

   Path target = to_dir_path / from_file_path.filename();
   fs::rename(from_file_path, target);
}

void
CopyFile(const Path& from_file_path, const Path& to_dir_path)
{
   ASSERT(FileExists(from_file_path), "The file ", from_file_path.filename(), " does not exist.")

   Path target = to_dir_path / from_file_path.filename();
   if(!FileExists(to_dir_path)) fs::create_directories(to_dir_path); // Recursively create target directory, if it doesn't exist.
   fs::copy_file(from_file_path, target, fs::copy_options::overwrite_existing);
}

/***************************************************************************************************************************************************************
* Directory Handling Functions
***************************************************************************************************************************************************************/
bool
isDirectory(const Path& path) { return fs::is_directory(path); }

bool
DirectoryExists(const Path& dir_path) { return FileExists(dir_path) && isDirectory(dir_path); }

void
ClearDirectory(const Path& dir_path) { FOR_EACH(entry, fs::directory_iterator(dir_path)) fs::remove_all(entry.path()); }

void
CreateDirectory(const Path& dir_path, const bool clear_if_exists)
{
   ASSERT(DirectoryExists(dir_path.parent_path()), "The parent directory of the given path ", dir_path.filename(), " does not exist.")

   fs::create_directory(dir_path);
   if(clear_if_exists) ClearDirectory(dir_path);
}

void
CreateDirectories(const Path& dir_path) { fs::create_directories(dir_path); }

bool
DirectoryIsEmpty(const Path& dir_path)
{
   ASSERT(DirectoryExists(dir_path), "The directory ", dir_path.filename(), " does not exist.")
   return FileIsEmpty(dir_path) && isDirectory(dir_path);
}

bool
DeleteDirectory(const Path& dir_path)
{
   ClearDirectory(dir_path);
   return DeleteFile(dir_path);
}

void
CopyDirectory(const Path& from_dir_path, const Path& to_dir_path)
{
   ASSERT(DirectoryExists(from_dir_path), "The directory ", from_dir_path.filename(), " does not exist.")

   Path target = to_dir_path / from_dir_path.filename();
   fs::create_directories(to_dir_path); // Recursively create target directory, if it doesn't exist.
   fs::copy(from_dir_path, target, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
}

void
MoveDirectory(const std::string_view& from_dir_path, const std::string_view& to_dir_path)
{
   CopyDirectory(from_dir_path, to_dir_path);
   DEBUG_WARN_IF(!DeleteDirectory(from_dir_path), "Failed to delete the directory ", Path(from_dir_path).filename(), " during move operation.")
}

/***************************************************************************************************************************************************************
* General Linux Commands
***************************************************************************************************************************************************************/
int
RunCommand(const std::string& cmd, const bool show_output)
{
   DEBUG_ASSERT(system(nullptr), "The command processor was not found.")
   return system((cmd + (!show_output ? " >> /dev/null" : "")).c_str());
}

int
RunCommandFrom(const std::string& cmd, const Path& dir_path, const bool show_output)
{
   ASSERT(isDirectory(dir_path), "The given path must be a directory.")
   return RunCommand("cd " + dir_path.string() + "; " + cmd, show_output);
}

void
RunCommands(const DArray<std::string>& cmds, const bool show_output) { FOR_EACH_CONST(cmd, cmds) RunCommand(cmd, show_output); }

void
RunCommandsFrom(const DArray<std::string>& cmds, const Path& dir_path, const bool show_output)
{
   FOR_EACH_CONST(cmd, cmds) RunCommandFrom(cmd, dir_path, show_output);
}

bool
CommandExists(const std::string& cmd)
{
   ASSERT(!isSubstring(" " , cmd), "Cannot have spaces in the Linux command: ", cmd)
   return RunCommand("which " + cmd) == 0;
}

void
CompileTeXFile(const std::string& tex_compiler, const Path& tex_path, bool show_output)
{
   ASSERT(CommandExists(tex_compiler), "Please install ", tex_compiler, " to compile a .tex file or use a different compiler.")
   ASSERT(tex_path.extension() == ".tex", "The given file has the extension ", tex_path.extension(), ". Expected a .tex file.")
   RunCommandFrom(tex_compiler + " " + tex_path.filename().string(), tex_path.parent_path(), show_output);
}

void
ConvertPDFtoPNG(Path pdf_path, const UInt pixel_density, bool show_output)
{
   ASSERT(CommandExists("convert"), "Please install Magick to convert a .pdf file to a .png file.")
   ASSERT(pdf_path.extension() == ".pdf", "The given file has the extension ", pdf_path.extension(), ". Expected a .pdf file.")

   std::string cmd;
   cmd += "convert -density " + ToString(pixel_density) + " ";
   cmd += pdf_path;
   cmd += " -channel RGBA ";
   cmd += " -colorspace sRGB ";
   cmd += "PNG32:";
   cmd += pdf_path.replace_extension(".png");

   RunCommand(cmd, show_output);
}

}