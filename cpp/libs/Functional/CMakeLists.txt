include_directories(${PROJECT_SOURCE_DIR}/libs/Functional)

set(SOURCE_FILES
        include/Explicit.h
        include/Piecewise.h
        src/Explicit.cpp
        src/Piecewise.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary
        LinearAlgebraLibrary)

add_library(FunctionalLibrary ${SOURCE_FILES})
target_link_libraries(FunctionalLibrary ${LINK_LIBRARIES})
