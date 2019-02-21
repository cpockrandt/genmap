# ===========================================================================
#  GenMap tests
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

include (CTest)

# include (ExternalProject)
# ExternalProject_Add ("googletest"
#                      GIT_REPOSITORY "https://github.com/google/googletest.git"
#                      GIT_TAG "release-1.8.1"
#                      CMAKE_ARGS "-DBUILD_GMOCK=0 -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/contrib"
#                      UPDATE_DISCONNECTED YES
# )

# add_executable(my_test tests.cpp)
#
# add_dependencies(my_test googletest)
# include_directories(${PROJECT_BINARY_DIR}/contrib/include)
# target_link_libraries(my_test ${PROJECT_BINARY_DIR}/contrib/lib/libgtest.a)
# target_link_libraries(my_test ${PROJECT_BINARY_DIR}/contrib/lib/libgtest_main.a)
# target_link_libraries(my_test pthread)

# add_test(NAME    my_test
#          COMMAND my_test)

# TODO: add dependencies

function(add_test_suite TEST_NAME_PREFIX TEST_CASE_FOLDER INDEX_FLAGS MAP_FLAGS)
    # add_test(
    #     NAME "${TEST_NAME_PREFIX}_raw_map"
    #     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r" "raw_map")
    # add_test(
    #     NAME "${TEST_NAME_PREFIX}_raw_freq8"
    #     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r -fs" "raw_freq8")
    # add_test(
    #     NAME "${TEST_NAME_PREFIX}_raw_freq16"
    #     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r -fl" "raw_freq16")
    add_test(
        NAME "${TEST_NAME_PREFIX}_txt_map"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -t" "txt_map")
    add_test(
        NAME "${TEST_NAME_PREFIX}_txt_freq16"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -t -fl" "txt_freq16")
    add_test(
        NAME "${TEST_NAME_PREFIX}_wig_map"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -w" "wig_map")
    add_test(
        NAME "${TEST_NAME_PREFIX}_wig_freq16"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -w -fl" "wig_freq16")
    add_test(
        NAME "${TEST_NAME_PREFIX}_bed_map"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -b" "bed_map")
    add_test(
        NAME "${TEST_NAME_PREFIX}_bed_freq16"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -b -fl" "bed_freq16")
    add_test(
        NAME "${TEST_NAME_PREFIX}_csv"
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -d" "csv")
endfunction(add_test_suite)

add_test_suite("single_fasta_single_sequence"                 "1a" "-F"  "-E 0 -K 3")
add_test_suite("single_fasta_single_sequence_rc"              "1b" "-F"  "-E 0 -K 3 -c")
add_test_suite("single_fasta_multi_sequence"                  "2a" "-F"  "-E 0 -K 4")
add_test_suite("single_fasta_multi_sequence_rc"               "2b" "-F"  "-E 0 -K 4 -c")
add_test_suite("single_fasta_multi_sequence_hard_raw"         "2c" "-F"  "-E 0 -K 4")
add_test_suite("single_fasta_multi_sequence_hard_raw_rc"      "2d" "-F"  "-E 0 -K 4 -c")
add_test_suite("multi_fasta_multi_sequence"                   "3a" "-FD" "-E 0 -K 4")
add_test_suite("multi_fasta_multi_sequence_rc"                "3b" "-FD" "-E 0 -K 4 -c")
add_test_suite("multi_fasta_multi_sequence_exclude_pseudo"    "3c" "-FD" "-E 0 -K 4 -ep")
add_test_suite("multi_fasta_multi_sequence_exclude_pseudo_rc" "3d" "-FD" "-E 0 -K 4 -ep -c")
