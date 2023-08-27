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
#include "GLTypes.h"

#include <functional>
#include <memory>
#include <optional>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace aprn::vis {

class Animator;

/***************************************************************************************************************************************************************
* Action Class Full Specialisations
***************************************************************************************************************************************************************/
template<ActionType T>
class Action : public ActionBase
{
 public:
   Action(Animator& animator)
      : ActionBase(animator, T, Zero, Zero, nullptr) { EXIT("The given action has not yet been implemented.") }

   void Do(const Real global_time) override {}
};

///** Ramp
//***************************************************************************************************************************************************************/
//template<ActionType T>
//requires Ramp<T>
//class Action<T> final : public ActionBase
//{
// public:
//   Action(Animator& animator, RampType _ramp_type, Real start_time, Real end_time, std::function<Real(Real)> reparam = Linear);
//
//   void Do(const Real global_time) override;
//
// private:
//   RampType Ramp{RampType::Scale};
//};

/** Offset
***************************************************************************************************************************************************************/
template<ActionType T>
requires Offset<T>
class Action<T> final : public ActionBase
{
 public:
   Action(Animator& animator, const glm::vec3& position_offset);

   Action(Animator& animator, Real angle_offset, const glm::vec3& axis);

   void Do(const Real global_time) override;

 private:
   glm::vec3 Position;
   Real      Angle;
   glm::vec3 Axis;
};

/** Scale
***************************************************************************************************************************************************************/
template<ActionType T>
requires Scale<T>
class Action<T> final : public ActionBase
{
 public:
   Action(Animator& animator, Real scale, Real start_time, Real end_time, Reparametriser reparam = Linear);

   Action(Animator& animator, const glm::vec3& scales, Real start_time, Real end_time, Reparametriser reparam = Linear);

   void Do(const Real global_time) override;

 private:
   const glm::vec3 Scales{One};
};

/** Translation
***************************************************************************************************************************************************************/
template<ActionType T>
requires Translation<T>
class Action<T> final : public ActionBase
{
 public:
   Action(Animator& animator, const glm::vec3& disp_or_posi, Real start_time, Real end_time, Reparametriser reparam = Linear);

   Action(Animator& animator, std::function<SVectorR3(Real)> path, Real start_time, Real end_time);

   void Do(const Real global_time) override;

 private:
   std::function<glm::vec3(Real)> Displacement;
};

/** Rotation
***************************************************************************************************************************************************************/
template<ActionType T>
requires Rotation<T>
class Action<T> final : public ActionBase
{
 public:
   Action(Animator& animator, Real angle, const glm::vec3& axis, Real start_time, Real end_time, Reparametriser reparam = Linear);

   Action(Animator& animator, Real angle, const glm::vec3& axis, const glm::vec3& ref_point, Real start_time, Real end_time,
          Reparametriser reparam = Linear);

   Action(Animator& animator, const glm::vec3& angular_velocity, Real start_time = Zero, Reparametriser ramp = Identity);

   void Do(const Real global_time) override;

 private:
   std::function<Real(Real)> Angle;
   const glm::vec3           Axis;
   const glm::vec3           Reference;
   Real                      AngularSpeed;
};

///** Reflection
//***************************************************************************************************************************************************************/
//template<ActionType T>
//requires Reflection<T>
//class Action<T> final : public ActionBase
//{
// public:
//   Action(Animator& animator, const SVectorF3& _normal, const SVectorF3& ref_point, Real start_time, Real end_time, std::function<Real(Real)> reparam = Linear);
//
//   void Do() override;
//
// private:
//   const SVectorF3 Normal;
//   const SVectorF3 ReferencePoint;
//   Action<ActionType::Scale> Scale;
//   Action<ActionType::MoveBy> Move;
//};
//

///** Morph
//***************************************************************************************************************************************************************/
//template<ActionType T>
//requires Morph<T>
//class Action<T> final : public ActionBase
//{
// public:
//   Action(Animator& animator, const Model& _other_actor, Real start_time, Real end_time, std::function<Real(Real)> reparam = Linear);
//
//   void Do() override;
//
// private:
//   const Model& OtherActor;
//   Action<ActionType::Scale> Scales;
//   Action<ActionType::MoveBy> Move;
//   Action<ActionType::RotateAbout> Rotate;
//};

///** Set Stroke/Fill Colour
//***************************************************************************************************************************************************************/
//template<ActionType T>
//requires SetColour<T>
//class Action<T> final : public ActionBase
//{
// public:
//   Action(Animator& animator, Real _new_colour, Real start_time, Real end_time, std::function<Real(Real)> reparam = Linear);
//
//   void
//   Do() override;
//
// private:
//   const SVectorF4 OldColour;
//   const SVectorF4 NewColour;
//};

///** Glow/Blink
//***************************************************************************************************************************************************************/
//template<ActionType T>
//requires SetColour<T>
//class Action<T> final : public ActionBase
//{
// public:
//   Action(Animator& animator, Real _new_colour, Real start_time, Real end_time, std::function<Real(Real)> reparam = Linear);
//
//   void
//   Do() override;
//
// private:
//   const SVectorF4 OldColour;
//   const SVectorF4 NewColour;
//};

}

#include "Action.tpp"
