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

#include "../include/Scene.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Scene Public Interface
***************************************************************************************************************************************************************/
Scene::Scene()
   : Scene(1000.0, true) {}

Scene::Scene(const Real duration, const bool adjust_duration)
   : Duration_(duration), AdjustDuration_(adjust_duration)
{
   ASSERT(Positive(duration) || adjust_duration, "Cannot have a negative duration for a scene unless the final duration is to be computed.")
   ASSERT(SingleScene_, "This constructor should only be called for the first scene.")

   SingleScene_ = false;
}

Scene::Scene(Scene& prev_scene, const Real duration, const bool adjust_duration)
   : Duration_(duration), AdjustDuration_(adjust_duration)
{
   ASSERT(Zero < duration      , "Cannot have a negative duration for a scene.")
   ASSERT(!SingleScene_        , "This constructor should not be called for the first scene.")
   ASSERT(PrevScene_           , "The current scene has already been assigned a previous scene.")
   ASSERT(prev_scene.NextScene_, "The previous scene has already been assigned a next scene.")

   PrevScene_ = &prev_scene;
   prev_scene.NextScene_ = this;
}

void
Scene::Init(const Real start_time)
{
   // Initialise all actorss except for those associated to TeXBoxes.
   FOR_EACH(actor, Actors_) if(!std::dynamic_pointer_cast<TeXBox>(actor)) actor->Init();

   // Initialise all lights.
   FOR_EACH(dlight, DLights_) dlight.Init();
   FOR_EACH(plight, PLights_) plight.Init();
   FOR_EACH(slight, SLights_) slight.Init();

   // Compute the start and end times of the scene.
   const auto max_duration(1.0e5);
   if(AdjustDuration_)
   {
      Real duration = -One;
      FOR_EACH_CONST(model, Actors_) if(model->ExitTime() < max_duration) Maximise(duration, model->ExitTime());
      Duration_ = Positive(duration) ? duration : Duration_;
      ASSERT(Positive(Duration_), "Could not adjust the scene duration based on model lifetimes. Please specify the duration for scene: ", Title_)
   }
   else
   {
      FOR_EACH_CONST(model, Actors_)
         if(model->ExitTime() < max_duration)
            ASSERT(model->ExitTime() < Duration_, "This model's lifespan exceeds that of scene: ", Title_)
   }

   StartTime_ = start_time;
   EndTime_   = StartTime_ + Duration_;
}

/***************************************************************************************************************************************************************
* Scene Private Interface
***************************************************************************************************************************************************************/
void
Scene::UpdateActors(const Real current_time) { FOR_EACH(actor, Actors_) actor->Update(current_time); }

void
Scene::RenderDirecShadows(Shader& shader)
{
   if(DLights_.empty()) return;
   else ASSERT(DLights_.size() == 1, "Can currently only handle one directional light.")

   shader.Bind();
   shader.SetDirectionalLightSpaceMatrix(DLights_[0].LightSpaceMatrix());

   auto& shadow_map = DLights_[0].ShadowMap();
   GLCall(glViewport(0, 0, shadow_map.DepthMap().Width(), shadow_map.DepthMap().Height()))

//  GLCall(glCullFace(GL_FRONT)); // Prevents peter-panning

   shadow_map.StartWrite();
   RenderModels(shader);
   shadow_map.StopWrite();

//   GLCall(glCullFace(GL_BACK));

   shader.Unbind();
}

void
Scene::RenderPointShadows(Shader& shader)
{
   if(PLights_.empty()) return;
   else ASSERT(PLights_.size() < 5, "Can currently handle at most four point lights.")

   shader.Bind();

   FOR_EACH(point_light, PLights_)
   {
      shader.SetPointLightSpaceMatrices(point_light.LightSpaceMatrices());
      shader.SetPointPosition(point_light.Position());
      shader.SetPointFarPlane(PointLight::FarPlane());

      auto& shadow_map = point_light.ShadowMap();
      GLCall(glViewport(0, 0, shadow_map.DepthMap().Width(), shadow_map.DepthMap().Height()))

      shadow_map.StartWrite();
      RenderModels(shader);
      shadow_map.StopWrite();
   }

   shader.Unbind();
}

void
Scene::RenderScene(Shader& shader, Camera& camera, const bool post_process)
{
   ASSERT(DLights_.size() <= 1, "Can currently only handle at most one directional light.")

   shader.Bind();

   shader.UseCamera(camera);
   shader.UsePostProcessor(post_process);

   if(!DLights_.empty())
   {
      shader.UseLight(DLights_[0]);
      shader.SetDirectionalLightSpaceMatrix(DLights_[0].LightSpaceMatrix());
      DLights_[0].ShadowMap().StartRead(1);
      shader.SetDirectionalShadowMap(1);
   }

   size_t i = 0;
   FOR_EACH(point_light, PLights_)
   {
      shader.UseLight(point_light);
      point_light.ShadowMap().StartRead(i + 2);
      shader.SetPointShadowMap(i, i + 2);
      i++;
   }
   shader.SetPointFarPlane(PointLight::FarPlane());

   RenderModels(shader);

   shader.Unbind();
}

void
Scene::RenderModels(Shader& shader)
{
   shader.SetUniform1i("u_use_diffuse_map"     , 0);
   shader.SetUniform1i("u_use_normal_map"      , 0);
   shader.SetUniform1i("u_use_displacement_map", 0);

   // Render model and its sub-models.
   FOR_EACH(actor, Actors_) actor->Render(shader);
}

}