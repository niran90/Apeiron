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
#include "ActionBase.h"
#include "Texture.h"

#include <functional>
#include <unordered_map>
#include <string>
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

class Shader;

class RenderObject
{
 public:
   virtual ~RenderObject() = default;

   virtual void Init() = 0;

   virtual bool Initialised() const = 0;

   virtual void ComputeLifespan() = 0;

   virtual void LoadTextureMap(const std::unordered_map<std::string, Texture&>& texture_map) = 0;

   virtual void Update(Real global_time) = 0;

   virtual void Render(Shader& shader) = 0;

   virtual void Delete() = 0;

   /** Set Object Attributes
   ************************************************************************************************************************************************************/
   virtual RenderObject& SetName(const std::string& name) = 0;

   virtual RenderObject& SetColour(const SVectorR4& rgba_colour) = 0;

   virtual RenderObject& SetColour(const Colour& colour) = 0;

   virtual RenderObject& SetMaterial(const std::string& name, Real specular_intensity, Real smoothness) = 0;

   virtual RenderObject& SetTexture(const std::string& material, size_t index, size_t resolution, Real dispacement_scale) = 0;

   virtual RenderObject& SetTexture(const std::string& material, const std::string& item, size_t index, size_t resolution, Real dispacement_scale) = 0;

   /** Set Object Actions
   ************************************************************************************************************************************************************/
   virtual RenderObject& OffsetPosition(const SVectorR3& displacement) = 0;

   virtual RenderObject& OffsetOrientation(Real angle, const SVectorR3& axis) = 0;

   virtual RenderObject& Scale(Real factor, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& Scale(const SVectorR3& factors, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& MoveBy(const SVectorR3& displacement, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& MoveTo(const SVectorR3& position, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& MoveAt(const SVectorR3& velocity, Real start_time, Reparametriser ramp) = 0;

   virtual RenderObject& Trace(std::function<SVectorR3(Real)> path, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& RotateBy(Real angle, const SVectorR3& axis, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& RotateAt(const SVectorR3& angular_velocity, Real start_time, Reparametriser ramp) = 0;

   virtual RenderObject& RevolveBy(Real angle, const SVectorR3& axis, const SVectorR3& refe_point, Real start_time, Real end_time, Reparametriser reparam) = 0;

   virtual RenderObject& RevolveAt(const SVectorR3& angular_velocity, const SVectorR3& refe_point, Real start_time, Reparametriser ramp) = 0;

   /** Accessors
   ************************************************************************************************************************************************************/
   inline const auto& Name() const { return Name_; }

   inline const auto& EntryTime() const { return EntryTime_; }

   inline const auto& ExitTime() const { return ExitTime_; }

 protected:
   std::string  Name_;
   Real         EntryTime_{Zero};
   Real         ExitTime_{InfFloat<>};
   mutable bool Init_{false};
};

}
