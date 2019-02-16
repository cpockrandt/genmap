# ===========================================================================
#  GenMap tests
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

include (ExternalProject)
ExternalProject_Add ("googletest"
                     GIT_REPOSITORY "https://github.com/google/googletest.git"
                     GIT_TAG "release-1.8.1"
                     CMAKE_ARGS "-DBUILD_GMOCK=0 -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/contrib"
                     UPDATE_DISCONNECTED YES
)

include (CTest)

add_executable(my_test tests.cpp)

add_dependencies(my_test googletest)
include_directories(${PROJECT_BINARY_DIR}/contrib/include)
target_link_libraries(my_test ${PROJECT_BINARY_DIR}/contrib/lib/libgtest.a)
target_link_libraries(my_test ${PROJECT_BINARY_DIR}/contrib/lib/libgtest_main.a)
target_link_libraries(my_test pthread)

add_test(NAME    my_test
         COMMAND my_test)

add_test(NAME test1
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh
         "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}")
add_dependencies(test1 genmap)
# foreach(DI sa fm)
#     add_test (NAME test_mkindex_${PROG}_${DI}
#               COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/maintests.sh
#                       "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" ${PROG} ${DI} "MKINDEX" " ")
# endforeach()
