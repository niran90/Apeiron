include_directories(${PROJECT_SOURCE_DIR}/libs/FileManager)

set(SOURCE_FILES
        include/File.h
        include/File.tpp
        include/FileSystem.h
        src/File.cpp
        src/FileSystem.cpp)

set(LINK_LIBRARIES)

add_library(FileManagerLibrary ${SOURCE_FILES})
target_link_libraries(FileManagerLibrary ${LINK_LIBRARIES})
