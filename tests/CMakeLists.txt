# ===========================================================================
#  GenMap tests
# ===========================================================================

include (ExternalProject)
ExternalProject_Add (googletest
                     PREFIX "${CMAKE_CURRENT_BINARY_DIR}/googletest"
                     GIT_REPOSITORY "https://github.com/google/googletest.git"
                     INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}"
                     CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
                     GIT_TAG release-1.8.0
                     UPDATE_DISCONNECTED YES)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") # TODO: --coverage
add_executable (genmap_algo_test tests.cpp)
add_dependencies (genmap_algo_test googletest)

include_directories (${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries (genmap_algo_test ${SEQAN_LIBRARIES})
target_link_libraries (genmap_algo_test ${CMAKE_CURRENT_BINARY_DIR}/lib/libgtest.a)
target_link_libraries (genmap_algo_test ${CMAKE_CURRENT_BINARY_DIR}/lib/libgtest_main.a)
target_link_libraries (genmap_algo_test pthread)

add_test(NAME algo_test COMMAND genmap_algo_test)

# TODO: add dependencies
function (add_test_suite TEST_NAME_PREFIX TEST_CASE_FOLDER INDEX_FLAGS MAP_FLAGS)
    add_test (NAME "${TEST_NAME_PREFIX}_raw_map"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r" "raw_map")
    add_test (NAME "${TEST_NAME_PREFIX}_raw_freq8"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r -fs" "raw_freq8")
    add_test (NAME "${TEST_NAME_PREFIX}_raw_freq16"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -r -fl" "raw_freq16")
    add_test (NAME "${TEST_NAME_PREFIX}_txt_map"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -t" "txt_map")
    add_test (NAME "${TEST_NAME_PREFIX}_txt_freq16"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -t -fl" "txt_freq16")
    add_test (NAME "${TEST_NAME_PREFIX}_wig_map"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -w" "wig_map")
    add_test (NAME "${TEST_NAME_PREFIX}_wig_freq16"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -w -fl" "wig_freq16")
    add_test (NAME "${TEST_NAME_PREFIX}_bed_map"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -b" "bed_map")
    add_test (NAME "${TEST_NAME_PREFIX}_bed_freq16"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -b -fl" "bed_freq16")
    add_test (NAME "${TEST_NAME_PREFIX}_csv"
              COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "${TEST_CASE_FOLDER}" "${INDEX_FLAGS}" "${MAP_FLAGS} -d" "csv")
endfunction (add_test_suite)

add_test_suite ("single_fasta_single_sequence"                 "1a" "-F"  "-E 0 -K 3")
add_test_suite ("single_fasta_single_sequence_rc"              "1b" "-F"  "-E 0 -K 3 -c")
add_test_suite ("single_fasta_multi_sequence"                  "2a" "-F"  "-E 0 -K 4")
add_test_suite ("single_fasta_multi_sequence_rc"               "2b" "-F"  "-E 0 -K 4 -c")
add_test_suite ("single_fasta_multi_sequence_hard_raw"         "2c" "-F"  "-E 0 -K 4")
add_test_suite ("single_fasta_multi_sequence_hard_raw_rc"      "2d" "-F"  "-E 0 -K 4 -c")
add_test_suite ("multi_fasta_multi_sequence"                   "3a" "-FD" "-E 0 -K 4")
add_test_suite ("multi_fasta_multi_sequence_rc"                "3b" "-FD" "-E 0 -K 4 -c")
add_test_suite ("multi_fasta_multi_sequence_exclude_pseudo"    "3c" "-FD" "-E 0 -K 4 -ep")
add_test_suite ("multi_fasta_multi_sequence_exclude_pseudo_rc" "3d" "-FD" "-E 0 -K 4 -ep -c")
