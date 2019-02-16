# ===========================================================================
#                  GenMap
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
string(ASCII 27 Esc)
set(ColourBold  "${Esc}[1m")
set(ColourReset "${Esc}[m")
set(ColourRed   "${Esc}[31m")

message ("${ColourBold}Compiler Detection${ColourReset}")

project (genmap CXX)

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Dependency detection${ColourReset}")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/include
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_PREFIX_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_PREFIX_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
    message (STATUS "Found a local SeqAn library provided with the GenMap source code.")
    message (    "   This will be preferred over system global headers.")
endif ()

# ----------------------------------------------------------------------------
# Add GenMap targets
# ----------------------------------------------------------------------------
add_subdirectory(src)

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Add Tests
# ----------------------------------------------------------------------------
message ("\n${ColourBold}Setting up unit tests${ColourReset}")
enable_testing()
add_subdirectory(tests)
