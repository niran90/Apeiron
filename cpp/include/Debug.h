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

#include <iostream>

namespace aprn{

/***************************************************************************************************************************************************************
* General Error Handling
***************************************************************************************************************************************************************/

/** Get a string with file name and line number. */
#define FILE_LINE(file, line) "\n FILE: " + ToString(file) + "\n LINE: " + ToString(line) + "\n"

/** Print a Ping from the current file and line. */
#define PING usleep(100); Print("Ping from file:", __FILE__, "at line number:", __LINE__, "\n"); Flush();

/** Throw error without exiting from a given file and line. */
#define ERROR_FROM(type, file, line, args...)\
{\
  std::cerr<<"\n=================================================================================================================",\
             "\n ", type, ": ", args, \
                 FILE_LINE(file, line), \
               "=================================================================================================================\n";\
}

/** Throw warning without exiting from a given file and line. */
#define WARNING_FROM(file, line, args...)\
{\
  std::cout<<"\n=================================================================================================================",\
             "\n WARN: ", args, \
                 FILE_LINE(file, line), \
               "=================================================================================================================\n";\
}

/** Throw warning (without exiting) from the current line and file. */
#define WARN(args...) WARNING_FROM(__FILE__, __LINE__, args)

/** Throw warning (without exiting) from the current line and file. */
#define WARN_IF(condition, args...) if(!(condition)) WARN(args)

/** Throw runtime error and exit. */
#define EXIT(args...) { ERROR_FROM("RUNTIME ERROR", __FILE__, __LINE__, args); exit(-1); }

/** Throw runtime error from a given file and line, and exit. */
#define EXIT_FROM(file, line, args...) { ERROR_FROM("RUNTIME ERROR", file, line, args); exit(-1); }

/** Throw compile-time error and exit. */
#define STATIC_EXIT(args...) { ERROR_FROM("COMPILE-TIME ERROR", __FILE__, __LINE__, args); exit(-1); }

/** Assert a static condition. If false, throw an error and exit. */
#define STATIC_ASSERT(static_condition, args...) if constexpr(!(static_condition)) STATIC_EXIT(args);

/** Assert a condition. If false, throw an error and exit. */
#define ASSERT(condition, args...) if(!(condition)) EXIT(args);

/** Debug mode warnings and assertions. */
#ifdef DEBUG_MODE
   #define DEBUG_WARN(args...)               WARN(args);
   #define DEBUG_WARN_IF(condition, args...) WARN_IF(condition, args);
   #define DEBUG_ASSERT(condition, args...)  ASSERT(condition, args);
#else
   #define DEBUG_WARN(args...)              {};
   #define DEBUG_WARN_IF(args...)           {};
   #define DEBUG_ASSERT(condition, args...) {};
#endif

/***************************************************************************************************************************************************************
* Numerical Error Handling
***************************************************************************************************************************************************************/

#define ASSERT_NUMBER(value) ASSERT(isNumber(value), "The passed value is not a number.");

#define DEBUG_ASSERT_NUMBER(value) DEBUG_ASSERT(isNumber(value), "The passed value is not a number.");

#define ASSERT_VALID_NUMBER(value) if(std::isnan(value) || std::isinf(value)) EXIT(std::fpclassify(ax) == 0 ? "NaN detected\n\t\t" : "Inf detected\n\t\t");

}