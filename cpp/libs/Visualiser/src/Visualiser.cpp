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

#include "../include/Visualiser.h"

#include <execution>
#include <optional>
#include <string>

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Public Interface
***************************************************************************************************************************************************************/
Visualiser::Visualiser()
//   : Visualiser(1920, 1080) {}
//   : Visualiser(2560, 1440) {}
   : Visualiser(3840, 2160) {}

Visualiser::Visualiser(const GLint window_width, const GLint window_height)
   : Window_(window_width, window_height), Cameras_{Camera()}, ActiveCamera_(&Cameras_.front()) {}

void
Visualiser::Add(Scene& scene) { Add(std::move(scene)); }

void
Visualiser::Add(const Camera& camera) { Cameras_.push_back(camera); }

void
Visualiser::Add(Scene&& scene) { Scenes_.push_back(std::move(scene)); }

void
Visualiser::Add(Camera&& camera) { Cameras_.push_back(std::move(camera)); }

void
Visualiser::Animate()
{
   // Initialise all visualiser elements.
   Init();

   // Run main render loop.
   while(Window_.isOpen())
   {
      BeginFrame();
      UpdateScene();
      HandleUserInputs();
      RenderScene();
      PostProcess();
      RenderGUIWindow();
      EndFrame();
   }

   // Destroy all graphics contexts.
   Terminate();
}

/***************************************************************************************************************************************************************
* Private Interface
***************************************************************************************************************************************************************/
void
Visualiser::Init()
{
   InitWindow();
   InitGUI();
   InitScenes();
   InitTeXBoxes();
   InitTextures();
   InitCameras();
   InitShaders();
   InitPostProcessor();
   Window_.InitTime();
}

void
Visualiser::InitWindow()
{
   // Open a window, set its title, and initialise OpenGL.
   Window_.Open();
   Window_.SetTitle("Apeiron Visualiser");
   Window_.InitOpenGL();
}

void
Visualiser::InitGUI()
{
#ifdef DEBUG_MODE
   GUI_.Init(Window_.GLFWWindow_);
#endif
}

void
Visualiser::InitScenes()
{
   auto first_scene_it = std::find_if(Scenes_.begin(), Scenes_.end(), [](const auto& scene){ return scene.PrevScene_ == nullptr; });
   ASSERT(first_scene_it != Scenes_.end(), "Failed to locate the first scene.")

   CurrentScene_ = &(*first_scene_it);
   Scene* current_scene(CurrentScene_);
   size_t scene_count{};
   Real   start_time{};

   // Loop through scenes
   do
   {
      // Update current scene if it isn't the first iteration
      if(scene_count)
      {
         // Sync the start-time of the next scene (or the transition to it) to the end-time of the current scene.
         start_time = current_scene->EndTime_;

         // Initialise the scene transition if there is one, and re-update the start time of the next scene.
         if(current_scene->Transition_._Type != TransitionType::None)
         {
            auto& transition = current_scene->Transition_;
            transition.Init(start_time);
            start_time = transition._EndTime;
         }
         current_scene = current_scene->NextScene_; // Move over to the next scene
      }

      // Initialise current scene and update count
      current_scene->Init(start_time);
      ++scene_count;
   }
   while(current_scene->NextScene_);

   ASSERT(scene_count == Scenes_.size(), "There was a mismatch in the total number of scenes.")
}

void
Visualiser::InitTeXBoxes()
{
   // Linearise pointers to all TeX-boxes to allow for parallel initialisation.
   DArray<TeXBox*> tex_boxes;
   tex_boxes.reserve(10 * Scenes_.size());
   FOR_EACH(scene, Scenes_) FOR_EACH(tex_box, scene.TeXBoxes_) tex_boxes.push_back(tex_box.get());

   // Initialise LaTeX compilation directory, compile all LaTeX source code, generate glyph sheets, and initialise underlying tex-box Model.
   InitTeXDirectory();
   FOR(i, tex_boxes.size()) tex_boxes[i]->InitTeXBox(i);

   // Load tex-box model textures. Note: only diffuse texture required.
   FOR_EACH(scene, Scenes_)
      FOR(i, scene.TeXBoxes_.size())
      {
         auto& tex_box = scene.TeXBoxes_[i];

         const auto texture_name = "TeXture_" + ToString(i);
         const auto texture_type = TextureType::Diffuse;
         const auto type_string  = TextureTypeString(texture_type);

         // Load compiled tex-box image as a texture.
         UMap<Texture> texture_files;
         texture_files.emplace(type_string, Texture(texture_type, tex_box->ImagePath()));
         Textures_.emplace(texture_name, std::move(texture_files));

         // Point to the texture from the tex-box sub-glyph models.
         UMap<Texture&> texture_file_map;
         FOR_EACH(sub_texture_name, sub_texture, Textures_[texture_name]) texture_file_map.emplace(sub_texture_name, sub_texture);
         tex_box->LoadTextureMap(texture_file_map);
      }
}

void
Visualiser::InitTextures()
{
   // Load model textures
   FOR_EACH(scene, Scenes_)
      FOR_EACH_CONST(actor, scene.Actors_)
         if(auto model = std::dynamic_pointer_cast<Model>(actor))
            if(model->TextureRequest())
            {
               const auto& texture_info = model->TextureRequest().value();
               const auto& texture_name = texture_info.first;

               if(!Textures_.contains(texture_name))
               {
                  // Add all files associated to the given texture
                  const auto texture_list = { TextureType::Diffuse,
                                              TextureType::Normal,
                                              TextureType::Displacement }; // Note: add appropriate enums if more textures are to be read.
                  UMap<Texture> texture_files;
                  FOR_EACH_CONST(texture_type, texture_list)
                  {
                     const auto path = TexturePath(TextureDirectory(texture_name), texture_type);
                     if(path)
                     {
                        texture_files.emplace(TextureTypeString(texture_type), Texture(texture_type, path.value()));
                        if(texture_type == TextureType::Displacement)
                        {
                           const auto& displacement_map_scale = texture_info.second;
                           texture_files.at(TextureTypeString(texture_type)).SetMapScale(displacement_map_scale);
                        }
                     }
                     else EXIT("Could not locate the texture files of texture ", texture_name)
                  }

                  // Add texture files to the list of textures
                  Textures_.emplace(texture_name, std::move(texture_files));
               }

               // Point to the textures from the model.
               UMap<Texture&> texture_file_map;
               FOR_EACH(sub_texture_name, sub_texture, Textures_[texture_name]) texture_file_map.emplace(sub_texture_name, sub_texture);
               model->LoadTextureMap(texture_file_map);
            }
}

void
Visualiser::InitCameras()
{
   // Set default settings of the current active camera.
   ASSERT(ActiveCamera_, "The active camera pointer has not yet been set.")
   ActiveCamera_->SetOrientation(glm::vec3(0.0f, 0.0f, 2.0f), 0.0, 0.0);
   ActiveCamera_->SetViewFrustum(Window_.ViewportAspectRatio(), 45.0, 1.0, -100.0);
}

void
Visualiser::InitShaders()
{
   for(std::string shader : {"Default", "DirecShadow", "PointShadow", "Line"})
      Shaders_.emplace(shader, Shader::Directory + shader + ".glsl");
}

void
Visualiser::InitPostProcessor() { if(PostProcess_) PostProcessor_.Init(Window_.ViewportDimensions()); }

void
Visualiser::BeginFrame()
{
   // Clear currently bound frame buffer.
   ClearFrameBuffer();

   // Set new GUI frame, if debugging.
#ifdef DEBUG_MODE
   GUI_.NewFrame();
#endif

   // Update the current and previous times, compute delta time, compute and display frame-rate, and check if the viewport was modified.
   Window_.ComputeDeltaTime();
   Window_.ComputeFrameRate();
   ViewPortModified_ = Window_.ViewportModified();
}

void
Visualiser::UpdateScene()
{
   DEBUG_ASSERT(CurrentScene_, "The current scene pointer has not yet been set.")

   // Check if the current scene needs to be updated to the next.
   const auto current_time = Window_.CurrentTime();
   if(!CurrentScene_->isCurrent(current_time))
   {
      ASSERT(CurrentScene_->NextScene_, "TODO - need to properly handle final scene.")
      CurrentScene_ = CurrentScene_->NextScene_;
   }

   // Update all actors in the current scene.
   CurrentScene_->UpdateActors(current_time);
}

void
Visualiser::HandleUserInputs()
{
   // Handle cursor, key, and mouse wheel inputs.
   if(HideCursor_) ActiveCamera_->CursorControl(Window_.CursorDisplacement());
   ActiveCamera_->KeyControl(Window_.Keys_, Window_.DeltaTime());
   ActiveCamera_->WheelControl(Window_.WheelDisplacement());

   // If the viewport was modified, update the view frustum and adjust line shader resolution.
   if(ViewPortModified_)
   {
      ActiveCamera_->SetViewFrustum(Window_.ViewportAspectRatio());
      Shaders_.at("Line").SetUniform2f("u_resolution", Window_.ViewportDimensions_[0], Window_.ViewportDimensions_[1]);
   }
}

void
Visualiser::RenderScene()
{
   // Render shadows from all directional and point light sources.
   CurrentScene_->RenderDirecShadows(Shaders_.at("DirecShadow"));
   CurrentScene_->RenderPointShadows(Shaders_.at("PointShadow"));

   // Point shadow rendering modifies the viewport, so need to reset it.
   Window_.ResetViewport();

   // Write to off-screen frame buffer.
   if(PostProcess_) PostProcessor_.StartWrite();

   // Render all elements of the current scene.
   CurrentScene_->RenderScene(Shaders_.at("Default"), *ActiveCamera_, PostProcess_);

   // Finalise off-screen render.
   if(PostProcess_) PostProcessor_.StopWrite();
}

void
Visualiser::PostProcess() { if(PostProcess_) PostProcessor_.Render(); }

void
Visualiser::RenderGUIWindow()
{
#ifdef DEBUG_MODE
   // Start a GUI window and add required elements to it.
   GUI_.StartWindow();
   AddGUIElements();
   GUI_.EndWindow();

   // Render the GUI elements in the window.
   GUI_.Render();
#endif
}

void
Visualiser::AddGUIElements()
{
#ifdef DEBUG_MODE
   using namespace ImGui;

   auto& scene = *CurrentScene_; // Only add current scene elements.

   if(CollapsingHeader("General"))
   {

   }

   if(CollapsingHeader("Models"))
   {

   }

   if(CollapsingHeader("TeXBoxes"))
   {

   }

   if(CollapsingHeader("Lights"))
   {
      if(TreeNode("Direct lights"))
      {
         TreePop();
         Separator();
      }

      if(TreeNode("Point lights"))
      {
         FOR_EACH(light, scene.PLights_)
            if(TreeNode(light.Name().c_str()))
            {
               light.AddGUIElements();
               TreePop();
               Separator();
            }
         TreePop();
         Separator();
      }

      if(TreeNode("Spot lights"))
      {
         TreePop();
         Separator();
      }
   }

   if(CollapsingHeader("Camera"))
   {

   }

   if(CollapsingHeader("Post-processor"))
   {

   }
#endif
}

void
Visualiser::EndFrame()
{
   Window_.SwapBuffers();
   glfwPollEvents();
}

void
Visualiser::Terminate()
{
#ifdef DEBUG_MODE
   GUI_.Terminate();
#endif
   Window_.Terminate();
}

}

