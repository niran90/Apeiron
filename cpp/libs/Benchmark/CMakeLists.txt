include_directories(${PROJECT_SOURCE_DIR}/libs/Benchmark)

add_library(BenchmarkLibrary INTERFACE)
target_include_directories(BenchmarkLibrary INTERFACE include/)
