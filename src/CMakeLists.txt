# ----------------------------------------------------------------------------
# Build configuration
# ----------------------------------------------------------------------------

message ("\n${ColourBold}Build configuration${ColourReset}")

# ----------------------------------------------------------------------------
# Make "Release" the default cmake build type
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# ----------------------------------------------------------------------------
# Warn if cmake build type is not "Release"
# ----------------------------------------------------------------------------

if (NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "${ColourRed}CMAKE_BUILD_TYPE is not \"Release\", your binaries will be slow.${ColourReset}")
endif ()

message (STATUS "GenMap version is: ${SEQAN_APP_VERSION}")

# option (LAMBDA_MMAPPED_DB     "Use mmapped access to the database."                                                    OFF)
option (GENMAP_NATIVE_BUILD   "Architecture-specific optimizations, i.e. g++ -march=native."                           ON)
# option (LAMBDA_MULTIOPT_BUILD "Build and install multiple binaries with different optimisation levels + a dispatcher." OFF)
option (GENMAP_STATIC_BUILD   "Include all libraries in the binaries."                                                 OFF)

if (GENMAP_NATIVE_BUILD)
    add_definitions (-DGENMAP_NATIVE_BUILD=1)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -xHOST -ipo -no-prec-div -fp-model fast=2")
    endif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
endif (GENMAP_NATIVE_BUILD)

if (GENMAP_STATIC_BUILD)
    add_definitions (-DGENMAP_STATIC_BUILD=1)
    set (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    # apple does not support fully static builds, but at least libgcc and libstdc++
    if (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        message (WARNING "WARNING: Builds on Mac are never fully static.")
    else (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    endif (APPLE)
    # on linux cmake adds -rdynamic automatically which clang can't handle in static builds
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
endif (GENMAP_STATIC_BUILD)

# if (LAMBDA_MMAPPED_DB)
#     add_definitions (-DLAMBDA_MMAPPED_DB=1)
# endif (LAMBDA_MMAPPED_DB)

message (STATUS "The following options are selected for the build:")
message (   "     GENMAP_NATIVE_BUILD             ${GENMAP_NATIVE_BUILD}")
message (   "     GENMAP_STATIC_BUILD             ${GENMAP_STATIC_BUILD}")
message (STATUS "Run 'cmake -LH' to get a comment on each option.")
message (STATUS "Remove CMakeCache.txt and re-run cmake with -DOPTIONNAME=ON|OFF to change an option.")

# strip binaries to make them smaller
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif ()

# TODO: verify whether this works for Clang!
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wl,-s")
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Enable global exception handler for all seqan apps.
set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1)

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

# Set the right output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic") # TODO: --coverage nur für debug auf travis!

# Update the list of file names below if you add source files to your application.
set (GENMAP_SOURCE_FILES genmap.cpp
                         common.hpp
                         indexing.hpp
                         mappability.hpp
                         algo.hpp
                         output.hpp)

add_executable (genmap ${GENMAP_SOURCE_FILES})
target_link_libraries (genmap ${SEQAN_LIBRARIES})

# if (LAMBDA_MULTIOPT_BUILD)
#     add_executable (lambda2-sse4 ${LAMBDA_SOURCE_FILES})
#     target_link_libraries (lambda2-sse4 ${SEQAN_LIBRARIES})
#     set_target_properties (lambda2-sse4 PROPERTIES COMPILE_FLAGS "-mmmx -msse -msse2 -msse3 -mssse3 -msse4 -mpopcnt")
# endif ()

# ----------------------------------------------------------------------------
# Man-pages
# ----------------------------------------------------------------------------

# Umbrella man-page
# add_custom_command (OUTPUT lambda2.1
#                     COMMAND lambda2 --export-help=man > lambda2.1
#                     DEPENDS lambda2)
# # searchn subcommand
# add_custom_command (OUTPUT lambda2-searchn.1
#                     COMMAND lambda2 searchn --export-help=man > lambda2-searchn.1
#                     DEPENDS lambda2)
# # searchp subcommand
# add_custom_command (OUTPUT lambda2-searchp.1
#                     COMMAND lambda2 searchp --export-help=man > lambda2-searchp.1
#                     DEPENDS lambda2)
# # mkindexn subcommand
# add_custom_command (OUTPUT lambda2-mkindexn.1
#                     COMMAND lambda2 mkindexn --export-help=man > lambda2-mkindexn.1
#                     DEPENDS lambda2)
# # mkindexp subcommand
# add_custom_command (OUTPUT lambda2-mkindexp.1
#                     COMMAND lambda2 mkindexp --export-help=man > lambda2-mkindexp.1
#                     DEPENDS lambda2)
#
# add_custom_target (manual
#                    ALL
#                    DEPENDS lambda2.1 lambda2-searchn.1 lambda2-searchp.1 lambda2-mkindexn.1 lambda2-mkindexp.1)

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# # Adapt to system paths
# include (GNUInstallDirs)
#
# # if (LAMBDA_MULTIOPT_BUILD)
#     # # Install lambda binaries into LIBEXECDIR
#     # install (TARGETS lambda2 lambda2-sse4
#     #          DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})
#     #
#     # # Install dispatcher script into BINDIR
#     # configure_file(../bin/lambda2.in
#     #                ${CMAKE_CURRENT_BINARY_DIR}/lambda2-dispatch
#     #                @ONLY)
#     # install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda2-dispatch
#     #          DESTINATION ${CMAKE_INSTALL_BINDIR}
#     #          RENAME lambda2
#     #          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
# # elseif ()
#     # Install lambda into BINDIR, usually ${PREFIX}/bin
#     install (TARGETS genmap
#              DESTINATION ${CMAKE_INSTALL_BINDIR})
# # endif ()
#
# # Install non-binary files for the package to DOCDIR, usually ${PREFIX}/share/doc/lambda2
# install (FILES ../LICENSE.rst
#                ../LICENSE-BSD.rst
#                ../LICENSE-AGPL3.rst
#                ../README.rst
#          DESTINATION ${CMAKE_INSTALL_DOCDIR})
#
# # Man pages into MANDIR, usually ${PREFIX}/share/man/man1 (or without share)
# # install (FILES ${CMAKE_CURRENT_BINARY_DIR}/lambda2.1
# #                ${CMAKE_CURRENT_BINARY_DIR}/lambda2-searchn.1
# #                ${CMAKE_CURRENT_BINARY_DIR}/lambda2-searchp.1
# #                ${CMAKE_CURRENT_BINARY_DIR}/lambda2-mkindexn.1
# #                ${CMAKE_CURRENT_BINARY_DIR}/lambda2-mkindexp.1
# #          DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# Information
# set (CPACK_PACKAGE_NAME "Genmap")
# set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "GenMap -- Fast and Exact Computation of Genome Mappability")
# set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.rst")
# set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.rst")
# set (CPACK_PACKAGE_VENDOR "Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>")
# set (CPACK_PACKAGE_CONTACT "${CPACK_PACKAGE_VENDOR}")
# set (CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_APP_VERSION_MAJOR}")
# set (CPACK_PACKAGE_VERSION_MINOR "${SEQAN_APP_VERSION_MINOR}")
# set (CPACK_PACKAGE_VERSION_PATCH "${SEQAN_APP_VERSION_PATCH}")
# set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION}")
# set (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")
#
# # Package format(s)
# if (CMAKE_SYSTEM_NAME MATCHES "Windows")
#     set(CPACK_GENERATOR "ZIP;NSIS")
# elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
#     set(CPACK_GENERATOR "ZIP;DragNDrop")
# elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
#     set(CPACK_GENERATOR "TBZ2")
# else()
#     set(CPACK_GENERATOR "TXZ")
# endif ()
#
# if (CMAKE_SYSTEM_NAME MATCHES "Linux")
#     set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB;RPM")
# endif ()
#
# # Package architecture
# if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
#     set(CMAKE_SYSTEM_PROCESSOR "x86_64")
#     set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
# endif ()
# if (CMAKE_CXX_FLAGS MATCHES "avx2")
#     set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_avx2")
# elseif (CMAKE_CXX_FLAGS MATCHES "sse4")
#     set (CMAKE_SYSTEM_PROCESSOR "${CMAKE_SYSTEM_PROCESSOR}_sse4")
# endif()
#
# # Include architecture in package name
# if (NOT DEFINED CPACK_SYSTEM_NAME)
#   set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
# endif (NOT DEFINED CPACK_SYSTEM_NAME)
#
# include (CPack)
