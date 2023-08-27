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
#include <DataContainer/include/Array.h>

namespace aprn{

template<typename T, unsigned N>
struct SortObject
{
   size_t Index;
   SArray<T, N> Values;

   SortObject() {}

   SortObject(const size_t index, const SArray<T, N>& values) : Index(index), Values(values) {}
};

template<typename T, unsigned N = 1>
class Sort
{
 public:
   Sort() = default;

   Sort(const size_t object_count) { Init(object_count); }

   inline void Init(const size_t object_count)
   {
      STATIC_ASSERT(isArithmetic<T>(), "Can only sort numerical data types currently.")
      SortObjects_.reserve(object_count);
   }

   inline void AddSortObject(const size_t index, const SArray<T, N>& values) { SortObjects_.emplace_back(index, values); }

   inline void SortAll(const bool _is_sort_in_ascending = true)
   {
      std::sort(SortObjects_.begin(), SortObjects_.end(), [](const SortObject<T, N>& a, const SortObject<T, N>& b)
                                                          {
                                                             FOR(i, N)
                                                             {
                                                                if     (a[i] < b[i]) return true;
                                                                else if(a[i] > b[i]) return false;
                                                                else continue;
                                                             }
                                                          });
      if(!_is_sort_in_ascending) std::reverse(SortObjects_.begin(), SortObjects_.end());
   }

   inline size_t GetIndex(const size_t i) const { return SortObjects_[i].Index_; }

   inline const SArray<T, N>& GetValues(const size_t i) const { return SortObjects_[i].Values; }

 private:
   DArray<SortObject<T, N>> SortObjects_;
};

}
