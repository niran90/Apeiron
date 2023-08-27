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

#include "../include/GUI.h"

namespace aprn::vis {

void
GUI::Init(GLFWwindow* window)
{
   // Setup ImGui binding
   IMGUI_CHECKVERSION();
   ImGui::CreateContext();
   ImGuiIO_ = &ImGui::GetIO();
   static_cast<void>(*ImGuiIO_);
   ImGui::StyleColorsDark();

   ImGui_ImplGlfw_InitForOpenGL(window, true);
   ImGui_ImplOpenGL3_Init("#version 460");

   isInitialised_ = true;
}

void
GUI::NewFrame()
{
   DEBUG_ASSERT(isInitialised_, "The GUI has not yet been initialised.")
   ImGui_ImplOpenGL3_NewFrame();
   ImGui_ImplGlfw_NewFrame();
   ImGui::NewFrame();
}

void
GUI::StartWindow()
{
   DEBUG_ASSERT(isInitialised_, "The GUI has not yet been initialised.")

   // Set window position and size.
   const ImGuiViewport* main_viewport = ImGui::GetMainViewport();
   ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 650, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
   ImGui::SetNextWindowSize(ImVec2(550, 680), ImGuiCond_FirstUseEver);

   ImGui::Begin("Visualiser Debugging Window");
}

void
GUI::EndWindow() { ImGui::End(); }

void
GUI::Render()
{
   DEBUG_ASSERT(isInitialised_, "The GUI has not yet been initialised.")
   ImGui::Render();
   ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void
GUI::Terminate()
{
   DEBUG_ASSERT(isInitialised_, "The GUI has not yet been initialised.")
   ImGui_ImplOpenGL3_Shutdown();
   ImGui_ImplGlfw_Shutdown();
   ImGui::DestroyContext();
}

}
