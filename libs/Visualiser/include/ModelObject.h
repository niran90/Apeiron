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

#include <functional>
#include <GL/glew.h>
#include <glm/glm.hpp>

namespace aprn::vis {

class ModelObject
{
 public:
   virtual ~ModelObject() = default;

   virtual void Update(Real global_time) = 0;

   virtual void Render() = 0;

   virtual void Delete() = 0;

   /** Set Object Attributes
   ************************************************************************************************************************************************************/
   virtual ModelObject& SetColour(const SVectorR4& rgba_colour) = 0;

   virtual ModelObject& SetColour(const Colour& colour) = 0;

   virtual ModelObject& SetMaterial(const std::string& name, Real specular_intensity, Real smoothness) = 0;

   virtual ModelObject& SetTexture(const std::string& material, size_t index, size_t resolution, Real dispacement_scale) = 0;

   virtual ModelObject& SetTexture(const std::string& material, const std::string& item, size_t index, size_t resolution, Real dispacement_scale) = 0;

   /** Set Object Actions
   ************************************************************************************************************************************************************/
   virtual ModelObject& OffsetPosition(const SVectorR3& displacement) = 0;

   virtual ModelObject& OffsetOrientation(Real angle, const SVectorR3& axis) = 0;

   virtual ModelObject& Scale(Real factor, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& Scale(const SVectorR3& factors, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& MoveBy(const SVectorR3& displacement, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& MoveTo(const SVectorR3& position, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& MoveAt(const SVectorR3& velocity, Real start_time = Zero, const std::function<Real(Real)>& ramp = Identity) = 0;

   virtual ModelObject& Trace(std::function<SVectorR3(Real)> path, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& RotateBy(Real angle, const SVectorR3& axis, Real start_time, Real end_time, const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& RotateAt(const SVectorR3& angular_velocity, Real start_time = Zero, const std::function<Real(Real)>& ramp = Identity) = 0;

   virtual ModelObject& RevolveBy(Real angle, const SVectorR3& axis, const SVectorR3& refe_point, Real start_time, Real end_time,
                                  const std::function<Real(Real)>& reparam) = 0;

   virtual ModelObject& RevolveAt(const SVectorR3& angular_velocity, const SVectorR3& refe_point, Real start_time = Zero,
                                  const std::function<Real(Real)>& ramp = Identity) = 0;

   /** Accessors
   ************************************************************************************************************************************************************/
   virtual bool isInitialised() const = 0;

   inline auto EntryTime() const { return EntryTime_; }

   inline auto ExitTime() const { return ExitTime_; }

 protected:
   friend class Scene;
   friend class ModelGroup;

   virtual void Init() = 0;

   virtual void ComputeLifespan() = 0;

   Real         EntryTime_{Zero};
   Real         ExitTime_{InfFloat<>};
   mutable bool isInitialised_{false};
};

}