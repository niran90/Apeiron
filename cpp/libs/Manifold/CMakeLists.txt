include_directories(${PROJECT_SOURCE_DIR}/libs/Manifold)

set(SOURCE_FILES
        include/Curve.h
        include/Curve.tpp
        include/Surface.h
        src/Surface.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary
        LinearAlgebraLibrary)

add_library(ManifoldLibrary ${SOURCE_FILES})
target_link_libraries(ManifoldLibrary ${LINK_LIBRARIES})
