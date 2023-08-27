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
#include "Functional/include/Explicit.h"
#include "Manifold/include/Curve.h"
#include "ActionBase.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <map>

namespace aprn::vis {

class Model;

class Animator
{
 public:
   Animator() = default;

   Animator(const Animator& other) = default;

   Animator(Animator&& other) noexcept;

   void Add(const ActionType action_type, SPtr<ActionBase> action_ptr);

   Pair<Real> ComputeLifespan();

   void Update(Real global_time);

   void Scale(const glm::vec3& factors);

   void Translate(const glm::vec3& displacement);

   void Rotate(const GLfloat angle, const glm::vec3& axis);

   /** Assignment Operators
   ************************************************************************************************************************************************************/
   Animator& operator=(const Animator& animator) = default;

   Animator& operator=(Animator&& animator) noexcept;

   /** Accessors
   ************************************************************************************************************************************************************/
   const glm::vec3 Displacement(const glm::vec3& position) const;

   inline const auto& ModelMatrix() const { return ModelMatrix_; }

 private:
   void Reset();

   using ActionMap = std::multimap<ActionType, SPtr<ActionBase>, ActionTypeComparator>;

   Model*    Model_;
   ActionMap Actions_;
   glm::mat4 ModelMatrix_{1.0f};
   glm::mat4 PastActions_{1.0f};
//   UPtr<mnfld::Curve<3>> Boundary;
};

}
