include_directories(${PROJECT_SOURCE_DIR}/libs/DataContainer)

add_library(DataContainerLibrary INTERFACE)
target_include_directories(DataContainerLibrary INTERFACE include/)