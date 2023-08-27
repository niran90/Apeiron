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

#include "../include/Animator.h"
#include "../include/Action.h"
#include "../include/Model.h"

namespace aprn::vis {

Animator::Animator(Animator&& other) noexcept { *this = std::move(other); }

void
Animator::Add(const ActionType action_type, SPtr<ActionBase> action_ptr) { Actions_.insert({action_type, action_ptr}); }

Pair<Real>
Animator::ComputeLifespan()
{
   Real entry_time = MaxFloat<>;
   Real exit_time  = LowestFloat<>;

   // Compute lifespan of all model actions.
   FOR_EACH(_, action, Actions_)
   {
      entry_time = Min(entry_time, action->StartTime_);
      exit_time  = Max(exit_time , action->EndTime_);
   }

   return { entry_time, exit_time };
}

void
Animator::Update(Real global_time)
{
   Reset();
   FOR_EACH(_, action, Actions_) action->Do(global_time);
}

void
Animator::Reset() { ModelMatrix_ = glm::mat4(1.0); }

void
Animator::Scale(const glm::vec3& factors) { ModelMatrix_ = glm::scale(ModelMatrix_, factors); }

void
Animator::Translate(const glm::vec3& displacement) { ModelMatrix_ = glm::translate(ModelMatrix_, displacement); }

void
Animator::Rotate(const GLfloat angle, const glm::vec3& axis) { ModelMatrix_ = glm::rotate(ModelMatrix_, angle, axis); }

Animator&
Animator::operator=(Animator&& animator) noexcept
{
   Model_       = std::move(animator.Model_);
   Actions_     = std::move(animator.Actions_);
   ModelMatrix_ = std::move(animator.ModelMatrix_);
   PastActions_ = std::move(animator.PastActions_);

   // Reset moved-from animator as it is now in an undefined state. Note: to avoid an infinite regress, we need to specifically invoke the copy assigment
   // operator here, NOT the move assignment operator.
   animator = Unmove(Animator());

   // Tricky: when moving actions, need to re-assign the 'Animator_' member of each action to the current animator.
   FOR_EACH(_, action, Actions_) action->Animator_ = this;

   return *this;
}

const glm::vec3
Animator::Displacement(const glm::vec3& position) const { return position - Model_->Centroid_; }

}
