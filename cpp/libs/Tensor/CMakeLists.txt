include_directories(${PROJECT_SOURCE_DIR}/libs/Tensor)

set(SOURCE_FILES
        include/Tensor.h
        src/Tensor.cpp)

set(LINK_LIBRARIES
        DataContainerLibrary
        LinearAlgebraLibrary)

add_library(TensorLibrary ${SOURCE_FILES})
target_link_libraries(TensorLibrary ${LINK_LIBRARIES})
