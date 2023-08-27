include_directories(${PROJECT_SOURCE_DIR}/libs/Polytope)

set(SOURCE_FILES
        include/Categories.h
        include/Polygon.h
        include/Polyhedron.h
        include/Polytope.h
        src/Polygon.cpp src/Categories.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary
        LinearAlgebraLibrary)

add_library(PolytopeLibrary ${SOURCE_FILES})
target_link_libraries(PolytopeLibrary ${LINK_LIBRARIES})
