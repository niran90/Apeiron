include_directories(${PROJECT_SOURCE_DIR}/libs/Graph)

set(SOURCE_FILES
        include/ADTree.h
        include/Graph.h
        include/KDTree.h
        include/Tree.h
        src/Graph.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary)

add_library(GraphLibrary ${SOURCE_FILES})
target_link_libraries(GraphLibrary ${LINK_LIBRARIES})
