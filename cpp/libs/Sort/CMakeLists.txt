include_directories(${PROJECT_SOURCE_DIR}/libs/DataContainer)

set(SOURCE_FILES
        include/Sort.h
        src/Sort.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary
        LinearAlgebraLibrary)

add_library(SortLibrary ${SOURCE_FILES})
target_link_libraries(SortLibrary ${LINK_LIBRARIES})
