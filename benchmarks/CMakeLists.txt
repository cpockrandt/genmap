# ===========================================================================
#  GenMap benchmarks
# ===========================================================================

include (ExternalProject)
ExternalProject_Add (googlebenchmark
                     PREFIX googlebenchmark
                     GIT_REPOSITORY "https://github.com/google/benchmark.git"
                     GIT_TAG "v1.4.1"
                     CMAKE_ARGS "-DBENCHMARK_ENABLE_TESTING=false;-DCMAKE_BUILD_TYPE=Release;-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                     UPDATE_DISCONNECTED yes)

add_executable (genmap_algo_benchmark bench.cpp)
add_dependencies (genmap_algo_benchmark googlebenchmark)
set_source_files_properties (bench.cpp PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_BINARY_DIR}/include")
target_link_libraries (genmap_algo_benchmark "-L${CMAKE_CURRENT_BINARY_DIR}/lib -lbenchmark")
target_link_libraries (genmap_algo_benchmark pthread) # required by google benchmark
