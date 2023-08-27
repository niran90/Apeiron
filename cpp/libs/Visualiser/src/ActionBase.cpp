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

#include "../include/ActionBase.h"
#include "../include/Model.h"

namespace aprn::vis {

bool
isTimeParametrised(ActionType type)
{
   return type == ActionType::RotateAt ||
          type == ActionType::MoveAt ||
          type == ActionType::Trace;
}

/***************************************************************************************************************************************************************
* Abstract Action Base Class
***************************************************************************************************************************************************************/
ActionBase::ActionBase(Animator& animator, const ActionType action_type)
   : ActionBase(animator, action_type, -One, -One, nullptr)
{
   ASSERT(action_type == ActionType::OffsetOrientation || action_type == ActionType::OffsetPosition, "This contructor only accepts static actions.")
}

ActionBase::ActionBase(Animator& animator, const ActionType action_type, const Real start_time, const Real end_time, const Reparametriser reparam)
   : Animator_(&animator), Type_(action_type), StartTime_(start_time), EndTime_(end_time),
     ParamNormaliser_(One / (EndTime_ - StartTime_)), Reparametriser_(reparam) {}

ActionBase::ActionBase(Animator& animator, const ActionType action_type, const Real start_time, const Reparametriser ramp)
   : Animator_(&animator), Type_(action_type), StartTime_(start_time), EndTime_(InfFloat<>), Ramp_(ramp) {}

Option<Real>
ActionBase::ComputeParameter(const Real global_time)
{
   const auto local_time = global_time - StartTime_;

   if(isTimeParametrised(Type_)) return std::make_optional(Negative(local_time) ? -One : local_time);
   else
   {
      const auto param =
         global_time <= EndTime_ ? std::make_optional(Negative(local_time) ? -One : Reparametriser_(ParamNormaliser_ * local_time)) : std::nullopt;

      if(param && Positive(param.value()))
      ASSERT((isBounded<true, true, true>(param.value(), Zero, One)), "The parameter must be in the [0, 1] range.")
      return param;
   }
}

}