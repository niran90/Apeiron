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

namespace aprn::vis {

enum class TransitionType
{
   None
};

class Scene;
class Transition
{
 public:
   Transition() = default;

   Transition(TransitionType type, Real duration)
      : _Type(type), _Duration(duration) {}

 private:
   friend class Visualiser;

   void Init(const Real start_time);

   TransitionType _Type{TransitionType::None};
   SPtr<Scene>    _FromScene;
   SPtr<Scene>    _ToScene;
   Real          _Duration{};
   Real          _StartTime;
   Real          _EndTime;
};

}