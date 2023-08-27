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
#include "../../DataContainer/include/Array.h"

#include <chrono>

namespace aprn {

/***************************************************************************************************************************************************************
* Time Units
***************************************************************************************************************************************************************/
enum class TimeUnit
{
   NanoSecond,
   MicroSecond,
   MilliSecond,
   Second
};

/***************************************************************************************************************************************************************
* Timer Class
***************************************************************************************************************************************************************/
class Timer
{
   friend class StopWatch;
   friend class Benchmark;

 public:
   /** Default constructor. */
   Timer() : TotalNanoSecondsLapsed(Zero), isRunning(false) {}

   /********** Class Methods **********/

   /** Reset timer. */
   inline void Reset()
   {
     isRunning = false;
     TotalNanoSecondsLapsed = Zero;
   }

   /** Start timer. */
   inline void Start()
   {
     isRunning = true;
     StartTime = std::chrono::high_resolution_clock::now();
   }

   /** Stop timer and accumulate total time. */
   inline void Stop()
   {
     StopTime = std::chrono::high_resolution_clock::now();
     TotalNanoSecondsLapsed += std::chrono::duration_cast<std::chrono::nanoseconds>(StopTime - StartTime).count();
     isRunning = false;
   }

   /** Get total time lapsed. */
   inline Real TotalLapTime(const TimeUnit& _time_units = TimeUnit::MilliSecond) const
   {
     switch(_time_units)
     {
       case TimeUnit::NanoSecond:  return TotalNanoSecondsLapsed;
       case TimeUnit::MicroSecond: return 1.0e-3 * TotalNanoSecondsLapsed;
       case TimeUnit::MilliSecond: return 1.0e-6 * TotalNanoSecondsLapsed;
       case TimeUnit::Second:      return 1.0e-9 * TotalNanoSecondsLapsed;
       default: EXIT("Time unit not recognised.")
     }
   }

 private:
   std::chrono::time_point<std::chrono::high_resolution_clock> StartTime;
   std::chrono::time_point<std::chrono::high_resolution_clock> StopTime;
   Real TotalNanoSecondsLapsed;
   bool isRunning;
};

/***************************************************************************************************************************************************************
* Stopwatch Class
***************************************************************************************************************************************************************/
class StopWatch : public Timer
{
 public:
   Real LapTimeMin;
   Real LapTimeMax;
   Real LapTimeMean;
   Real LapTimeRMS;
   Real LapTimeStd;
   bool isFinalised;
   DynamicArray<Real> LapTimes;

   /********** Constructors/Destructors **********/

   /** Default constructor. */
   StopWatch()
      : LapTimeMin(MaxReal<>), LapTimeMax(LowestReal<>), LapTimeMean(Zero), LapTimeRMS(Zero), LapTimeStd(Zero), isFinalised(false)
   {
      LapTimes.reserve(1e5);
   }

   /** Default destructor. */
   ~StopWatch() = default;

   /********** Class Methods **********/

   /** Accumulate min, max, and lap times. */
   inline void AccumulateLapTimes(const TimeUnit& _time_units = TimeUnit::MilliSecond)
   {
      ASSERT(!isRunning, "Cannot accumulate - the timer is still running.")

      const Real total_time_lapsed = TotalLapTime(_time_units);
      LapTimeMin = Min(LapTimeMin, total_time_lapsed);
      LapTimeMax = Max(LapTimeMax, total_time_lapsed);
      LapTimes.push_back(total_time_lapsed);
   }

   /** Finalise the time metrics for the current lap. */
   inline void FinaliseLap()
   {
      ASSERT(!isRunning, "Cannot finalise - the timer is still running.")

      if(isFinalised) return;

      // Finalise the mean and RMS.
      LapTimeMean = std::accumulate(LapTimes.begin(), LapTimes.end(), Zero); // Note: Last argument determines the type of the return value.
      LapTimeMean /= (Real)LapTimes.size();

      // Finalise the RMS and standard deviation.
      LapTimeRMS = Zero;
      LapTimeStd = Zero;
      FOR_EACH(lap_time, LapTimes)
      {
        LapTimeRMS += iPow(lap_time, 2);
        LapTimeStd += iPow(lap_time - LapTimeMean, 2);
      }
      LapTimeRMS = Sqrt( LapTimeRMS / static_cast<Real>(LapTimes.size()) );
      LapTimeStd = Sqrt( LapTimeStd / static_cast<Real>(LapTimes.size()) );

      isFinalised = true;
   }

   /** Reset stopwatch to start a new lap. */
   inline void Reset()
   {
      Timer::Reset();
      isFinalised = false;
   }
};

}
