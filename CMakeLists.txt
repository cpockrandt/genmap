# ===========================================================================
#                  GenMap
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string (ASCII 27 Esc)
set (ColourBold  "${Esc}[1m")
set (ColourReset "${Esc}[m")
set (ColourRed   "${Esc}[31m")

project (genmap CXX)

# ----------------------------------------------------------------------------
# App version
# ----------------------------------------------------------------------------

# change this after every release
set (SEQAN_APP_VERSION_MAJOR "0")
set (SEQAN_APP_VERSION_MINOR "9")
set (SEQAN_APP_VERSION_PATCH "0")
# don't change the following
set (SEQAN_APP_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")
# adapt when necessary
set (MINIMUM_SEQAN_VERSION "2.4.0")

# ----------------------------------------------------------------------------
# Compiler detection
# ----------------------------------------------------------------------------

message ("${ColourBold}Compiler Detection${ColourReset}")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.1)
        message (FATAL_ERROR "Your GCC version is too old. Minimum version is GCC-4.9.1!")
        return ()
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
        message (FATAL_ERROR "Your Clang version is too old. Please upgrade to 3.8.0 or use GCC.")
        return ()
    elseif (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.8)
        message (WARNING "Your Clang version is too old, you will not have parallelism! Upgrade to 3.8.0 or newer.")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17)
        message (WARNING "Your Intel Compiler version is too old. Please upgrade to 17.0.0 or newer!")
    endif ()
else ()
    message (WARNING "Unknown compiler, you are own your own!")
endif ()

message (STATUS "Compiling with ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}.")

message (STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Compiler specifics
# ----------------------------------------------------------------------------

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ftemplate-depth-1024")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # parallelize parts of build even for one translation unit
    if (NOT DEFINED GENMAP_COMPILE_THREADS)
        include (ProcessorCount)
        ProcessorCount (GENMAP_COMPILE_THREADS)
    endif ()
    if (GENMAP_COMPILE_THREADS GREATER 1)
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto=${GENMAP_COMPILE_THREADS}")
    endif()

    # required, otherwise #pragma clang diagnostics will trigger warnings
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
else ()
    # TODO: parallelize for Clang (isn't this done by default by the default clang linker?)
    set (dummy_ignore "${GENMAP_COMPILE_THREADS}") # surpress warnings
endif ()

if (NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # -Wextra not met by current SeqAn on intel
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wextra")
endif ()

# ----------------------------------------------------------------------------
# Dependency detection
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Dependency detection${ColourReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
                            ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
                            ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
                            ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the GenMap source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

# Search SeqAn and select dependencies.
find_package (OpenMP QUIET)
find_package (SeqAn  QUIET REQUIRED CONFIG)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

message (STATUS "These dependencies where found:")
message (   "     OPENMP     ${OPENMP_FOUND}      ${OpenMP_CXX_FLAGS}")
message (   "     SEQAN      ${SEQAN_FOUND}      ${SEQAN_VERSION_STRING}")

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
    message (WARNING "WARNING WARNING WARNING\nWARNING: OpenMP not found. GenMap will be built without multi-threading! "
    "This is probably not what you want! Use GCC >= 4.9.1, Clang >= 3.8.0 or ICC >= 16.0.2\nWARNING WARNING WARNING")
endif (NOT OPENMP_FOUND)

if (SEQAN_VERSION_STRING VERSION_LESS "${MINIMUM_SEQAN_VERSION}")
    message (FATAL_ERROR "The minimum SeqAn verison required is ${MINIMUM_SEQAN_VERSION}!")
    return ()
endif ()

# deactivate the version check on broken seqan releases
if (SEQAN_VERSION_STRING VERSION_LESS "2.3.2")
    add_definitions (-DSEQAN_DISABLE_VERSION_CHECK="YES")
endif ()

message (STATUS "The requirements were met.")

# ----------------------------------------------------------------------------
# Add targets
# ----------------------------------------------------------------------------

add_subdirectory (src)

# ----------------------------------------------------------------------------
# Add tests
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Setting up unit tests${ColourReset}")
enable_testing ()
add_subdirectory (tests EXCLUDE_FROM_ALL)

message (STATUS "Tests added.")

# ----------------------------------------------------------------------------
# Add benchmarks
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Setting up benchmarks${ColourReset}")
add_subdirectory (benchmarks EXCLUDE_FROM_ALL)

message (STATUS "Benchmarks added.")

message ("\n${ColourBold}Finishing up${ColourReset}")
