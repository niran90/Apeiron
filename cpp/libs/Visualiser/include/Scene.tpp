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

template<class T>
Scene&
Scene::Add(T&& object)
{
   using type = RemoveConstRef<T>;

   // Note: move both lvalue and rvalue objects.
   if constexpr(isTypeSame<type, Model>() || isTypeSame<type, ModelGroup>()) Actors_.push_back(std::make_shared<type>(std::move(object)));
   else if constexpr(isTypeSame<type, DirectLight>()) DLights_.push_back(std::move(object));
   else if constexpr(isTypeSame<type, PointLight>())  PLights_.push_back(std::move(object));
   else if constexpr(isTypeSame<type, SpotLight>())   SLights_.push_back(std::move(object));
   else if constexpr(isTypeSame<type, TeXBox>())
   {
      TeXBoxes_.push_back(std::make_shared<type>(std::move(object)));
      Actors_.push_back(TeXBoxes_.back());
   }
   else EXIT("Currently cannot add the following object type: ", typeid(type).name())

   return *this;
}

}