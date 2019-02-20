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

set (CASE_1A_FLAGS "-E 0 -K 3")
# add_test(
#     NAME single_fasta_single_sequence_raw_map
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -r" "raw_map")
add_test(
    NAME single_fasta_single_sequence_raw_freq8
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -r -fs" "raw_freq8")
add_test(
    NAME single_fasta_single_sequence_raw_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_single_sequence_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -t" "txt_map")
add_test(
    NAME single_fasta_single_sequence_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_single_sequence_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -w" "wig_map")
add_test(
    NAME single_fasta_single_sequence_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_single_sequence_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -b" "bed_map")
add_test(
    NAME single_fasta_single_sequence_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -b -fl" "bed_freq16")
add_test(
    NAME single_fasta_single_sequence_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1a" "${CASE_1A_FLAGS} -d" "csv")

set (CASE_1B_FLAGS "-E 0 -K 3 -c")
# add_test(
#     NAME single_fasta_single_sequence_rc_raw_map
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -r" "raw_map")
# add_test(
#     NAME single_fasta_single_sequence_rc_raw_freq8
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -r -fs" "raw_freq8")
# add_test(
#     NAME single_fasta_single_sequence_rc_raw_freq16
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -t" "txt_map")
add_test(
    NAME single_fasta_single_sequence_rc_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -w" "wig_map")
add_test(
    NAME single_fasta_single_sequence_rc_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -b" "bed_map")
add_test(
    NAME single_fasta_single_sequence_rc_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -b -fl" "bed_freq16")
add_test(
    NAME single_fasta_single_sequence_rc_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "1b" "${CASE_1B_FLAGS} -d" "csv")

set (CASE_2A_FLAGS "-E 0 -K 4")
# add_test(
#     NAME single_fasta_multi_sequence_raw_map
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -r" "raw_map")
# add_test(
#     NAME single_fasta_multi_sequence_raw_freq8
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -r -fs" "raw_freq8")
# add_test(
#     NAME single_fasta_multi_sequence_raw_freq16
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_multi_sequence_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -t" "txt_map")
add_test(
    NAME single_fasta_multi_sequence_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_multi_sequence_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -w" "wig_map")
add_test(
    NAME single_fasta_multi_sequence_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_multi_sequence_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -b" "bed_map")
add_test(
    NAME single_fasta_multi_sequence_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -b -fl" "bed_freq16")
add_test(
    NAME single_fasta_multi_sequence_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2a" "${CASE_2A_FLAGS} -d" "csv")

set (CASE_2B_FLAGS "-E 0 -K 4 -c")
# add_test(
#     NAME single_fasta_multi_sequence_rc_raw_map
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -r" "raw_map")
# add_test(
#     NAME single_fasta_multi_sequence_rc_raw_freq8
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -r -fs" "raw_freq8")
# add_test(
#     NAME single_fasta_multi_sequence_rc_raw_freq16
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_multi_sequence_rc_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -t" "txt_map")
add_test(
    NAME single_fasta_multi_sequence_rc_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_multi_sequence_rc_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -w" "wig_map")
add_test(
    NAME single_fasta_multi_sequence_rc_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_multi_sequence_rc_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -b" "bed_map")
add_test(
    NAME single_fasta_multi_sequence_rc_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -b -fl" "bed_freq16")
add_test(
    NAME single_fasta_multi_sequence_rc_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2b" "${CASE_2B_FLAGS} -d" "csv")

set (CASE_2C_FLAGS "-E 0 -K 4")
# add_test(
#     NAME single_fasta_multi_sequence_hard_raw_map
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -r" "raw_map")
# add_test(
#     NAME single_fasta_multi_sequence_hard_raw_freq8
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -r -fs" "raw_freq8")
# add_test(
#     NAME single_fasta_multi_sequence_hard_raw_freq16
#     COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -r -fl" "raw_freq16")
add_test(
    NAME single_fasta_multi_sequence_hard_txt_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -t" "txt_map")
add_test(
    NAME single_fasta_multi_sequence_hard_txt_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -t -fl" "txt_freq16")
add_test(
    NAME single_fasta_multi_sequence_hard_wig_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -w" "wig_map")
add_test(
    NAME single_fasta_multi_sequence_hard_wig_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -w -fl" "wig_freq16")
add_test(
    NAME single_fasta_multi_sequence_hard_bed_map
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -b" "bed_map")
add_test(
    NAME single_fasta_multi_sequence_hard_bed_freq16
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -b -fl" "bed_freq16")
add_test(
    NAME single_fasta_multi_sequence_hard_csv
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests.sh "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" "2c" "${CASE_2C_FLAGS} -d" "csv")
