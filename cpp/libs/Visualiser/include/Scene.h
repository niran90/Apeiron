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
#include "TeXGlyph.h"
#include "Light.h"
#include "RenderObject.h"
#include "Model.h"
#include "ModelGroup.h"
#include "Shader.h"
#include "TeXBox.h"
#include "Transition.h"

#include <map>
#include <unordered_map>

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Scene Class
***************************************************************************************************************************************************************/
class Scene
{
 public:
   Scene();

   Scene(Real duration, bool adjust_duration = true);

   Scene(Scene& prev_scene, Real duration = -One, bool adjust_duration = true);

   Scene(const Scene& other) = delete;

   Scene(Scene&& other) noexcept = default;

   template<class T>
   Scene& Add(T&& object);

   inline bool isCurrent(const Real current_time) const { return StartTime_ <= current_time && current_time < EndTime_; }

 private:
   friend class Visualiser;
   friend class Transition;

   void Init(Real start_time);

   void UpdateActors(const Real current_time);

   void RenderDirecShadows(Shader& shader);

   void RenderPointShadows(Shader& shader);

   void RenderScene(Shader& shader, Camera& camera, bool post_process);

   void RenderModels(Shader& shader);

   template<class T> using UMap = std::unordered_map<std::string, T>;

   std::string                Title_;
   DArray<SPtr<RenderObject>> Actors_;
   DArray<SPtr<TeXBox>>       TeXBoxes_;
   DArray<DirectLight>        DLights_;
   DArray<PointLight>         PLights_;
   DArray<SpotLight>          SLights_;
   UMap<UMap<Texture&>>       Textures_;
   Transition                 Transition_{};
   Scene*                     PrevScene_{};
   Scene*                     NextScene_{};
   Real                       Duration_;
   Real                       StartTime_;
   Real                       EndTime_;
   bool                       AdjustDuration_{false};
   inline static bool         SingleScene_{true};
};

}

#include "Scene.tpp"