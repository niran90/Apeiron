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

#include "Animator.h"

namespace aprn::vis {

/***************************************************************************************************************************************************************
* Action Class Full Specialisations
***************************************************************************************************************************************************************/

///** Ramp
//***************************************************************************************************************************************************************/
//template<ActionType type>
//requires Ramp<type>
//Action<type>::Action(Animator& animator, RampType _ramp_type, Real start_time, Real end_time, std::function<Real(Real)> reparam)
//   : ActionBase(model, type, start_time, end_time, reparam), Ramp(_ramp_type) {}
//
//template<ActionType type>
//requires Ramp<type>
//void Action<type>::Do(const Real global_time)
//{
//
//}

/** Offset
***************************************************************************************************************************************************************/
template<ActionType type>
requires Offset<type>
Action<type>::Action(Animator& animator, const glm::vec3& position_offset)
   : ActionBase(animator, type), Position(position_offset) {}

template<ActionType type>
requires Offset<type>
Action<type>::Action(Animator& animator, Real angle_offset, const glm::vec3& axis)
   : ActionBase(animator, type), Angle(angle_offset), Axis(axis) {}

template<ActionType type>
requires Offset<type>
void
Action<type>::Do(const Real global_time)
{
   if     (type == ActionType::OffsetPosition)    Animator_->Translate(Position);
   else if(type == ActionType::OffsetOrientation) Animator_->Rotate(Angle, Axis);
}

/** Scale
***************************************************************************************************************************************************************/
template<ActionType type>
requires Scale<type>
Action<type>::Action(Animator& animator, Real scale, Real start_time, Real end_time, Reparametriser reparam)
   : ActionBase(animator, type, start_time, end_time, reparam), Scales(scale) {}

template<ActionType type>
requires Scale<type>
Action<type>::Action(Animator& animator, const glm::vec3& scales, Real start_time, Real end_time, Reparametriser reparam)
   : ActionBase(animator, type, start_time, end_time, reparam), Scales(scales) {}

template<ActionType type>
requires Scale<type>
void
Action<type>::Do(const Real global_time)
{
   const auto param = this->ComputeParameter(global_time);
   if(param && Positive(param.value())) Animator_->Scale(static_cast<float>(param.value()) * Scales);
   else {} // Remove action
}

/** Translation
***************************************************************************************************************************************************************/
template<ActionType type>
requires Translation<type>
Action<type>::Action(Animator& animator, const glm::vec3& disp_or_posi, Real start_time, Real end_time, Reparametriser reparam)
   : ActionBase(animator, type, start_time, end_time, reparam)
{
   switch(type)
   {
      case ActionType::MoveBy: this->Displacement = [&disp_or_posi](Real t){ return static_cast<float>(t) * disp_or_posi; }; break;
      case ActionType::MoveTo: this->Displacement = [&animator, &disp_or_posi](Real t){ return static_cast<float>(t) * animator.Displacement(disp_or_posi); }; break;
      default: throw std::invalid_argument("Unrecognised action type.");
   }
}

template<ActionType type>
requires Translation<type>
Action<type>::Action(Animator& animator, std::function<SVectorR3(Real)> path, Real start_time, Real end_time)
   : ActionBase(animator, type, start_time, end_time, [](Real){return Zero;})
{
   if(type == ActionType::Trace) this->Displacement = [path](Real t){ return SVectorToGlmVec(path(t)); };
   else throw std::invalid_argument("Unrecognised action type.");
}

template<ActionType type>
requires Translation<type>
void
Action<type>::Do(const Real global_time)
{
   const auto param = this->ComputeParameter(global_time);
   if(param && Positive(param.value())) Animator_->Translate(this->Displacement(param.value()));
   else {} // Remove action
}

/** Rotation
***************************************************************************************************************************************************************/
template<ActionType type>
requires Rotation<type>
Action<type>::Action(Animator& animator, Real angle, const glm::vec3& axis, Real start_time, Real end_time, Reparametriser reparam)
   : Action<type>::Action(animator, angle, axis, {0.0, 0.0, 0.0}, start_time, end_time, reparam) {}

template<ActionType type>
requires Rotation<type>
Action<type>::Action(Animator& animator, Real angle, const glm::vec3& axis, const glm::vec3& ref_point, Real start_time, Real end_time,
                     Reparametriser reparam)
   : ActionBase(animator, type, start_time, end_time, reparam), Axis(axis), Reference(ref_point)
{
   this->Angle = [angle](Real t){ return t * angle; };
}

template<ActionType type>
requires Rotation<type>
Action<type>::Action(Animator& animator, const glm::vec3& angular_velocity, Real start_time, Reparametriser ramp)
   : ActionBase(animator, type, start_time, ramp), Axis(glm::normalize(angular_velocity)), Reference({0.0, 0.0, 0.0}),
     AngularSpeed(glm::length(angular_velocity))
{
   this->Angle = [this](Real t){ return t * AngularSpeed * Ramp_(t); };
}

template<ActionType type>
requires Rotation<type>
void
Action<type>::Do(const Real global_time)
{
   const auto param = this->ComputeParameter(global_time);
   if(param && Positive(param.value()))
   {
      const bool is_revolution = type == ActionType::RevolveBy || type == ActionType::RevolveAt;
      if(is_revolution) Animator_->Translate(-Reference);
      Animator_->Rotate(this->Angle(param.value()), Axis);
      if(is_revolution) Animator_->Translate(Reference);
   }
   else {} // Remove action
}

///** Reflection
//***************************************************************************************************************************************************************/
//template<ActionType type>
//requires Reflect<type>
//Action<type>::Action(Animator& animator, const SVectorF3& _normal, const SVectorF3& ref_point, Real start_time, Real end_time,
//                     Reparametriser reparam)
//   : ActionBase(animator, type, start_time, end_time, reparam), Normal(_normal), ReferencePoint(ref_point)
//{
//
//}
//
//template<ActionType type>
//requires Reflect<type>
//void Action<type>::Do()
//{
//
//}
//

///** Morph
//***************************************************************************************************************************************************************/
//template<ActionType type>
//requires Morph<type>
//Action<type>::Action(Animator& animator, const Model& _other_actor, Real start_time, Real end_time, std::function<Real()> global_time,
//                     Reparametriser reparam)
//   : ActionBase(animator, start_time, end_time, global_time, reparam)
//{
//
//}
//
//template<ActionType type>
//requires Morph<type>
//void Action<type>::Do()
//{
//
//}
//
///** Set Stroke/Fill Colour
//***************************************************************************************************************************************************************/
//template<ActionType type>
//requires SetColour<type>
//Action<type>::Action(Animator& animator, Real _new_colour, Real start_time, Real end_time, std::function<Real()> global_time,
//                     Reparametriser reparam)
//   : ActionBase(animator, start_time, end_time, global_time, reparam)
//{
//
//}
//
//template<ActionType type>
//requires SetColour<type>
//void Action<type>::Do()
//{
//
//}

}

