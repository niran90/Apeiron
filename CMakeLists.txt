#***************************************************************************************************************************************************************
# Configure compile settings
#***************************************************************************************************************************************************************
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(Apeiron LANGUAGES C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_BUILD_TYPE Debug) # [Debug, Release]

# Set C++ compiler standard
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

# Set default, debug, and release flags
set(CMAKE_CXX_FLAGS         "-fopenmp -Wall -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -Wall -Wextra -Werror -pedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option")

#***************************************************************************************************************************************************************
# Setup build, library, and bin directories
#***************************************************************************************************************************************************************
set(BUILD_DIRECTORY "${CMAKE_SOURCE_DIR}/build")
set(LIB_DIRECTORY   "${BUILD_DIRECTORY}/lib")
set(BIN_DIRECTORY   "${BUILD_DIRECTORY}/bin")

# Set all compiler outputs to the build directory.
set(DEFAULT_OUT_OF_SOURCE_FOLDER ${BUILD_DIRECTORY})
set(CMAKE_BUILD_FILES_DIRECTORY  ${BUILD_DIRECTORY})
set(CMAKE_BUILD_DIRECTORY        ${BUILD_DIRECTORY})

# Set cmake library and binary output folders
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_PATH      ${LIB_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIRECTORY})

#***************************************************************************************************************************************************************
# Build instructions for OpenGL and its dependencies.
#***************************************************************************************************************************************************************

# Find and link OpenGL dependencies.
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OpenGL_FOUND)
    message("-- Error: OpenGL not found")
endif(NOT OpenGL_FOUND)

# Find and link GLEW.
find_package(GLEW 2.1.0 REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_directories(${GLEW_LIBRARY_DIRS})
add_definitions(${GLEW_DEFINITIONS})
if(NOT GLEW_FOUND)
    message("-- Error: GLEW not found")
endif(NOT GLEW_FOUND)

# Find and link GLFW.
find_package(glfw3 REQUIRED)
include_directories(${glfw3_INCLUDE_DIRS})
link_directories(${glfw3_LIBRARY_DIRS})
add_definitions(${glfw3_DEFINITIONS})
if(NOT glfw3_FOUND)
    message("-- Error: GLFW not found")
endif(NOT glfw3_FOUND)
set(GLFW_BUILD_DOCS OFF CACHE BOOL     "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL    "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

# Find and link GLM.
find_package(glm REQUIRED)
include_directories(${glm_INCLUDE_DIRS})
link_directories(${glm_LIBRARY_DIRS})
add_definitions(${glm_DEFINITIONS})
if(NOT glm_FOUND)
    message("Error: GLM not found")
endif(NOT glm_FOUND)

#***************************************************************************************************************************************************************
# Build instructions for ImGui and its dependencies.
#***************************************************************************************************************************************************************

# Compile ImGui as static library
set(IMGUI_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external/imgui)
file(GLOB IMGUI_SOURCE_FILES ${IMGUI_INCLUDE_DIRS}/*.cpp ${IMGUI_INCLUDE_DIRS}/backends/*.cpp)
add_library(ImGui STATIC ${IMGUI_SOURCE_FILES})
target_include_directories(ImGui SYSTEM PUBLIC ${IMGUI_INCLUDE_DIRS})
target_compile_options(ImGui PRIVATE -Wno-old-style-cast -Wno-sign-promo -Wno-strict-overflow -Wno-cast-qual)
target_link_libraries(ImGui ${CMAKE_DL_LIBS} glfw)

#***************************************************************************************************************************************************************
# Build instructions for OpenMP.
#***************************************************************************************************************************************************************
find_package(OpenMP REQUIRED)
include_directories(${OpenMP_INCLUDE_DIRS})
link_directories(${OpenMP_LIBRARY_DIRS})
add_definitions(${OpenMP_DEFINITIONS})
if(NOT OpenMP_FOUND)
    message("Error: OpenMP not found")
endif(NOT OpenMP_FOUND)

#***************************************************************************************************************************************************************
# Setup project sub-directories.
#***************************************************************************************************************************************************************
include_directories(${PROJECT_SOURCE_DIR}/libs)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/DataContainer)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Benchmark)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/FileManager)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Functional)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Graph)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/LinearAlgebra)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Manifold)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Polytope)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Sort)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Tensor)
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/Visualiser)

#***************************************************************************************************************************************************************
# Apeiron binary build instructions.
#***************************************************************************************************************************************************************
set(APEIRON_SOURCE_FILES
        src/main.cpp)

set(APEIRON_LINK_LIBRARIES
        DataContainerLibrary
        BenchmarkLibrary
        FileManagerLibrary
        FunctionalLibrary
        GraphLibrary
        LinearAlgebraLibrary
        ManifoldLibrary
        PolytopeLibrary
        SortLibrary
        TensorLibrary
        VisualiserLibrary)

add_executable(${PROJECT_NAME} ${APEIRON_SOURCE_FILES})

target_link_libraries(${PROJECT_NAME}
        ${APEIRON_LINK_LIBRARIES}
        ${OPENMP_LIBRARIES}
        ${OPENGL_LIBRARIES}
        ${GLEW_LIBRARIES}
        ${GLFW_LIBRARIES}
        ${GLM_LIBRARIES}
        ImGui)

#***************************************************************************************************************************************************************
# GoogleTest build instructions.
#***************************************************************************************************************************************************************
enable_testing()

# Fetch GoogleTest directories
include(FetchContent)
FetchContent_Declare(GoogleTest URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip)
FetchContent_MakeAvailable(GoogleTest)

include(CTest)
include(GoogleTest)

# Add unit test executables
add_executable(UnitTestBasicMath        ${PROJECT_SOURCE_DIR}/test/UnitTestBasicMath.cpp)
add_executable(UnitTestComparators      ${PROJECT_SOURCE_DIR}/test/UnitTestComparators.cpp)
add_executable(UnitTestMath             ${PROJECT_SOURCE_DIR}/test/UnitTestMath.cpp)
add_executable(UnitTestString           ${PROJECT_SOURCE_DIR}/test/UnitTestString.cpp)

add_executable(UnitTestArray            ${PROJECT_SOURCE_DIR}/libs/DataContainer/test/UnitTestArray.cpp)
add_executable(UnitTestNumericContainer ${PROJECT_SOURCE_DIR}/libs/DataContainer/test/UnitTestNumericContainer.cpp)
add_executable(UnitTestFileHandler      ${PROJECT_SOURCE_DIR}/libs/FileManager/test/UnitTestFileHandler.cpp)
add_executable(UnitTestParseTeX         ${PROJECT_SOURCE_DIR}/libs/Visualiser/test/UnitTestParseTeX.cpp)
add_executable(UnitTestVector           ${PROJECT_SOURCE_DIR}/libs/LinearAlgebra/test/UnitTestVector.cpp)
add_executable(UnitTestCurve            ${PROJECT_SOURCE_DIR}/libs/Manifold/test/UnitTestCurve.cpp)

# Link with gtest, gtest_main, and associated libraries.
target_link_libraries(UnitTestBasicMath        gtest gtest_main)
target_link_libraries(UnitTestComparators      gtest gtest_main)
target_link_libraries(UnitTestMath             gtest gtest_main)
target_link_libraries(UnitTestString           gtest gtest_main)

target_link_libraries(UnitTestArray            gtest gtest_main DataContainerLibrary)
target_link_libraries(UnitTestNumericContainer gtest gtest_main DataContainerLibrary)
target_link_libraries(UnitTestFileHandler      gtest gtest_main FileManagerLibrary)
target_link_libraries(UnitTestVector           gtest gtest_main LinearAlgebraLibrary)
target_link_libraries(UnitTestCurve            gtest gtest_main ManifoldLibrary)
target_link_libraries(UnitTestParseTeX         gtest gtest_main VisualiserLibrary)

# Add tests with CTest
gtest_discover_tests(UnitTestBasicMath)
gtest_discover_tests(UnitTestComparators)
gtest_discover_tests(UnitTestMath)
gtest_discover_tests(UnitTestString)
gtest_discover_tests(UnitTestArray)
gtest_discover_tests(UnitTestNumericContainer)
gtest_discover_tests(UnitTestFileHandler)
gtest_discover_tests(UnitTestVector)
gtest_discover_tests(UnitTestCurve)
gtest_discover_tests(UnitTestParseTeX)
