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

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Public Interface
***************************************************************************************************************************************************************/
Visualiser::Visualiser()
   : Visualiser(1920, 1080) {}

Visualiser::Visualiser(GLint window_width, GLint window_height)
   : _OpenGLWindow(window_width, window_height), _Cameras{{"Main", Camera()}}, _ActiveCamera(&_Cameras["Main"]) {}

void
Visualiser::Add(Scene& scene, const std::string& name)
{
   const std::string& id = name.empty() ? "Scene_" + ToStr(_Scenes.size()) : name;
   _Scenes.emplace(id, std::move(scene));
}

void
Visualiser::Add(Camera&& camera, const std::string& name)
{
   const std::string& id = name.empty() ? "Camera_" + ToStr(_Cameras.size()) : name;
   _Cameras.emplace(id, std::move(camera));
}

void
Visualiser::Render()
{
   Init();



   FrameBuffer fbo;
   RenderBuffer rbo;
   Texture cbo(TextureType::Diffuse, true);

   cbo.Init(1920, 1080, GL_RGB, GL_RGB, GL_UNSIGNED_BYTE, GL_CLAMP_TO_BORDER);
   rbo.Init();
   rbo.Allocate(GL_DEPTH24_STENCIL8, 1920, 1080);

   fbo.Init();
   fbo.Bind();
   fbo.AttachTexture2D(GL_COLOR_ATTACHMENT0, cbo.ID());
   fbo.AttachRenderBuffer(GL_DEPTH_STENCIL_ATTACHMENT, rbo.ID());
   fbo.Unbind();



   while(_OpenGLWindow.isOpen())
   {
      StartFrame();
      UpdateScene();
      ManageUserInputs();
      UpdateViewFrustum();
      RenderScene();
      EndFrame();
   }
}

/***************************************************************************************************************************************************************
* Private Interface
***************************************************************************************************************************************************************/
void
Visualiser::Init()
{
   // Set default settings of the main camera
   ASSERT(_ActiveCamera, "The active camera pointer has not yet been set.")
   _ActiveCamera->SetOrientation(glm::vec3(0.0f, 0.0f, 1.0f), 0.0, 90.0);
   _ActiveCamera->SetViewFrustum(_OpenGLWindow.ViewportAspectRatio(), 45.0, 1.0, -100.0);

   // Load all shaders
   _Shaders.emplace("General", "libs/Visualiser/resources/shaders/General.glsl");
   _Shaders.emplace("Line", "libs/Visualiser/resources/shaders/Line.glsl");
   _Shaders.emplace("DirectionalShadow", "libs/Visualiser/resources/shaders/DirectionalShadow.glsl");
   _Shaders.emplace("PointShadow", "libs/Visualiser/resources/shaders/PointShadow.glsl");

   // Initialise scenes, tex-boxes, and textures.
   InitScenes();
   InitTeXBoxes();
   InitTextures();

   // Set window title and set clock time to zero
   _OpenGLWindow.SetTitle("Apeiron");
   _OpenGLWindow.ResetTime();
}

void
Visualiser::InitScenes()
{
   auto first_scene_it = std::find_if(_Scenes.begin(), _Scenes.end(), [](const auto& entry){ return entry.second._PrevScene == nullptr; });
   ASSERT(first_scene_it != _Scenes.end(), "Could not locate the first scene.")

   _CurrentScene = &first_scene_it->second;
   Scene* current_scene(_CurrentScene);
   size_t scene_count{};
   Float  start_time{};

   // Loop through scenes
   do
   {
      // Update current scene if it isn't the first iteration
      if(scene_count)
      {
         // Sync the start-time of the next scene (or the transition to it) to the end-time of the current scene.
         start_time = current_scene->_EndTime;

         // Initialise the scene transition if there is one, and re-update the start time of the next scene.
         if(current_scene->_Transition._Type != TransitionType::None)
         {
            auto& transition = current_scene->_Transition;
            transition.Init(start_time);
            start_time = transition._EndTime;
         }
         current_scene = current_scene->_NextScene; // Move over to the next scene
      }

      // Initialise current scene and update count
      current_scene->Init(start_time);
      ++scene_count;
   }
   while(current_scene->_NextScene);

   ASSERT(scene_count == _Scenes.size(), "There was a mismatch in the total number of scenes.")
}

void
Visualiser::InitTeXBoxes()
{
   // Linearise pointers to all TeX-boxes to allow for parallel initialisation.
   DArray<std::pair<size_t, TeXBox*>> tex_boxes;
   tex_boxes.reserve(10 * _Scenes.size());
   FOR_EACH(_, scene, _Scenes) FOR_EACH(_, tex_box, scene._TeXBoxes) tex_boxes.push_back({ tex_boxes.size(), tex_box.get() });

   // Initialise LaTeX compilation directory, compile all LaTeX source code, and generate glyph sheets.
   TeXBox::InitTeXDirectory();
   FOR(i, tex_boxes.size()) tex_boxes[i].second->Init(i);

   // Load model textures. Note: only diffuse texture required.
   FOR_EACH(_, scene, _Scenes)
      FOR_EACH_CONST(_, tex_box, scene._TeXBoxes)
      {
         const auto texture_name = tex_box->_Label + "_texture";
         const auto texture_type = TextureType::Diffuse;

         UMap<Texture> texture_files;
         texture_files.emplace(TextureTypeString(texture_type), Texture(texture_type, tex_box->ImagePath()));
         _Textures.emplace(texture_name, std::move(texture_files));

         // Point to the textures from the scene
         UMap<Texture&> texture_file_map;
         FOR_EACH(sub_texture_name, sub_texture, _Textures[texture_name]) texture_file_map.emplace(sub_texture_name, sub_texture);
         scene._Textures.emplace(texture_name, texture_file_map);
      }
}

void
Visualiser::InitTextures()
{
   // Load model textures
   FOR_EACH(_, scene, _Scenes)
      FOR_EACH_CONST(_, model, scene._Models)
         if(model->_Texture.has_value() && !_Textures.contains(model->_Texture.value()))
         {
            // Add all files associated to the given texture
            const auto& texture_name = model->_Texture.value();
            const auto  texture_list = { TextureType::Diffuse,
                                         TextureType::Normal,
                                         TextureType::Displacement }; // Add appropriate enums if more textures are to be read
            UMap<Texture> texture_files;
            FOR_EACH_CONST(texture_type, texture_list)
            {
               const auto path = TexturePath(TextureDirectory(texture_name), texture_type);
               if(path.has_value())
               {
                  texture_files.emplace(TextureTypeString(texture_type), Texture(texture_type, path.value()));
                  if(texture_type == TextureType::Displacement)
                  {
                     const Float displacement_map_scale = 0.08; // TODO: Temporarily hard-code
                     texture_files.at(TextureTypeString(texture_type)).SetMapScale(displacement_map_scale);
                  }
               }
               else EXIT("Could not locate the texture files of texture ", texture_name)
            }

            // Add texture files to the list of textures
            _Textures.emplace(texture_name, std::move(texture_files));

            // Point to the textures from the scene
            UMap<Texture&> texture_file_map;
            FOR_EACH(sub_texture_name, sub_texture, _Textures[texture_name]) texture_file_map.emplace(sub_texture_name, sub_texture);
            scene._Textures.emplace(texture_name, texture_file_map);
         }
}

void
Visualiser::StartFrame()
{
   // Clear window
   GLCall(glClearColor(0.0f, 0.0f, 0.0f, 1.0f));
   GLCall(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));

   // Update the current and previous times, compute delta time, compute and display frame-rate, and check if the viewport was modified.
   _OpenGLWindow.ComputeDeltaTime();
   _OpenGLWindow.ComputeFrameRate();
   _isViewPortModified = _OpenGLWindow.isViewportModified();
}

void
Visualiser::EndFrame()
{
   _OpenGLWindow.SwapBuffers();
   glfwPollEvents();
}

void
Visualiser::UpdateScene()
{
   DEBUG_ASSERT(_CurrentScene, "The current scene pointer has not yet been set.")

   // Determine if the current scene needs to be updated
   const auto current_time = _OpenGLWindow.CurrentTime();
   if(!_CurrentScene->isCurrent(current_time)) _CurrentScene = _CurrentScene->_NextScene;

   // Update models in current scene
   _CurrentScene->UpdateModels(current_time);
}

void
Visualiser::ManageUserInputs()
{
   _ActiveCamera->KeyControl(_OpenGLWindow._Keys, _OpenGLWindow.DeltaTime());
   _ActiveCamera->CursorControl(_OpenGLWindow.CursorDisplacement());
   _ActiveCamera->WheelControl(_OpenGLWindow.WheelDisplacement());
}

void
Visualiser::UpdateViewFrustum()
{
   // If the viewport was modified, update the view frustum and adjust line shader resolution
   if(_isViewPortModified)
   {
      _ActiveCamera->SetViewFrustum(_OpenGLWindow.ViewportAspectRatio());
      _Shaders["Line"].SetUniform2f("u_resolution", _OpenGLWindow._ViewportDimensions[0], _OpenGLWindow._ViewportDimensions[1]);
   }
}

void
Visualiser::RenderScene()
{
   _CurrentScene->RenderDirectionalShadows(_Shaders["DirectionalShadow"]);
   _CurrentScene->RenderPointShadows(_Shaders["PointShadow"]);
   _OpenGLWindow.ResetViewport();
   _CurrentScene->RenderScene(_Shaders["General"], *_ActiveCamera);
}

}

