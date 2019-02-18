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

# single_fasta_single_seq
# single_fasta_multiple_seq
# directory_single_fasta_single_seq
# directory_multiple_fasta_mixed_seq

# TODO: add dependencies
add_test(
    NAME single_fasta_single_sequence_raw_freq8
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-r -fs" "raw_freq8")
add_test(
    NAME single_fasta_single_sequence_raw_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-r -fl" "raw_freq16")
add_test(
    NAME single_fasta_single_sequence_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-t" "txt_map")
add_test(
    NAME single_fasta_single_sequence_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-t -fl" "txt_freq16")
add_test(
    NAME single_fasta_single_sequence_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-w" "wig_map")
add_test(
    NAME single_fasta_single_sequence_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-w -fl" "wig_freq16")
add_test(
    NAME single_fasta_single_sequence_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-b" "bed_map")
add_test(
    NAME single_fasta_single_sequence_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-b -fl" "bed_freq16")
add_test(
    NAME single_fasta_single_sequence_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1" "-d" "csv")

# add_test(
#     NAME single_fasta_single_sequence_rc_raw_freq8
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -r -fs" "raw_freq8")
# add_test(
#     NAME single_fasta_single_sequence_rc_raw_freq16
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -t" "txt_map")
add_test(
    NAME single_fasta_single_sequence_rc_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -w" "wig_map")
add_test(
    NAME single_fasta_single_sequence_rc_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -b" "bed_map")
add_test(
    NAME single_fasta_single_sequence_rc_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -b -fl" "bed_freq16")
# add_test(
#     NAME single_fasta_single_sequence_rc_csv
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "-c -d" "csv")
